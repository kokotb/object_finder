import numpy as np
import cv2
import skimage.morphology
from scipy import ndimage
import skimage.measure
import time
import specpy as sp
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image
import os

from skimage.measure import label
from skimage.measure import regionprops



def get_initial_z(meas, config_num):
    """
    Gets the current z position
    """
    zrow_start, len_data, z_start = get_z(meas, config_num)
    if zrow_start == len_data:
        print('no signal!! stop the measurement!!')

    return zrow_start, z_start


def get_nucl_cent(nucl_m, img_part):
    all_centers = []
    nucl_labeled = label(nucl_m)
    dim = np.array(np.shape(nucl_m))
    min_img = dim * img_part
    max_img = (1 - img_part) * dim

    properties = regionprops(nucl_labeled)
    centers_x = []
    centers_y = []
    sceleton = []
    circ = []
    for props in properties:
        y0, x0 = props.centroid
        center_arr = [y0, x0]
        area = props.area
        perimeter = props.perimeter
        circularity = (4 * np.pi * area) / (perimeter ** 2)
        circ.append(circularity)
        if circularity > 0.7 and area > 100 and all(center_arr > min_img) and all(center_arr < max_img):
            centers_x.append(x0)
            centers_y.append(y0)
            circ.append(circularity)
            all_centers.append(center_arr)

    centers_fin = np.array(all_centers)
    np.random.shuffle(centers_fin)
    # plt.imshow(nucl_labeled)
    #
    # for x, y, t in zip(centers_x, centers_y, circ):
    #     plt.text(int(x), int(y), str(int(t*100)), fontsize=7, color='red')
    # plt.axis('off')
    # plt.show()

    return centers_fin


def get_z(meas, config_num):
    """
    Function calculates the z position of glass under the sample.
    It calculets the biggest difference in intensity of a few consecutive lines.

    :param meas: object of imspector measurement
    :param config_num: consecutive number of selected configuration
    """
    conf = meas.configuration(config_num)        # chose the configuration that has the xz (or yz) scan, for focus
    meas.activate(conf)
    im.run(meas)
    stack1 = conf.stack(0)              # chose the stack of configuration (which channel)
    data_stack = stack1.data()
    pos = meas.parameters('ExpControl/scan/range/coarse_z/g_off')

    # select all the parameters
    num_of_rows = 2             # number of rows you add
    num_of_stacks = 2000          # how many times you repeat the addition (how many stacks)
    mult_thresh = 2            # multiplication of threshold, where does the signal begin, needs adjusting each time
    num_row = 0                         # leave alone
    all_dif = []                        # leave alone
    sel_row = 0
    # get the threshold sum
    threshold_row = np.sum(data_stack[0, 0, 2])                 # sum of firs row, noise
    if threshold_row < 3:
        threshold_row = 3

    # go over whole stack and find where the signal begins
    for row in data_stack[0, 0]:
        row_sum1 = np.sum(row)
        if row_sum1 > threshold_row * mult_thresh:              # needs adjusting!
            sel_row = num_row
            break
        num_row += 1

    if num_row == len(data_stack[0, 0]):                    # if there is no signal higher than the trashold
        print('no signal!!!')
        chosen_row = num_row
        return chosen_row, len(data_stack[0, 0]), pos

    else:
        # sel_row = 0                       #if you want to start at the top of the image
        stack_sum1 = 0
        for i in range(sel_row, sel_row + num_of_rows):
            stack_sum1 = stack_sum1 + np.sum(data_stack[0, 0, i])
        stack_sum2 = 0
        for i in range(sel_row + num_of_rows, sel_row + 2 * num_of_rows):
            stack_sum2 = stack_sum1 + np.sum(data_stack[0, 0, i])
        diff = stack_sum2 - stack_sum1
        cur_row = sel_row + 2 * num_of_rows
        chosen_row = sel_row + 2 * num_of_rows

        # get all the stacks
        for stac in range(num_of_stacks):
            if (cur_row + num_of_rows) < len(data_stack[0, 0]):      # so that you do not go out range
                stack_sum1 = stack_sum2
                stack_sum2 = 0
                for j in range(cur_row, cur_row + num_of_rows):
                    stack_sum2 = stack_sum2 + np.sum(data_stack[0, 0, j])
                diff2 = stack_sum2 - stack_sum1
                all_dif.append(diff2)
                if diff2 > diff:
                    diff = diff2
                    chosen_row = j
                cur_row = cur_row + num_of_rows
            else:
                break

    return chosen_row, len(data_stack[0, 0]), pos


def adjust_focus(meas, config_num, z_start, zrow_start, max_z):
    """
    Function adjusts the z position of objective according to the new z position

    :param meas: object of imspector measurement
    :param config_num: consecutive number of selected configuration
    :param z_start: value of referential z position - z of objective (in m)
    :param zrow_start: value of the referential z position - which row in image
    :param max_z: maximum value of z (in m)
    """
    zrow_new, len_data, z_new = get_z(meas, config_num)
    conf = meas.configuration(config_num)
    psz = conf.parameters('ExpControl/scan/range/z/psz')
    len_z = conf.parameters('ExpControl/scan/range/z/res') * psz
    if zrow_new == len_data:
        print('no signal!!!')
        meas.set_parameters('ExpControl/scan/range/coarse_z/g_off',  z_start)
        time.sleep(1)
    else:
        delt_z_pix = zrow_new - zrow_start
        delt_z = delt_z_pix * psz
        new_z = z_new + delt_z
        time.sleep(1)

        if (new_z - z_start) < (
                200e-6 - len_z / 2) and new_z < max_z:      # condition for z (working distance for 60wi = 280e-6m)
            meas.set_parameters('ExpControl/scan/range/coarse_z/g_off', new_z)  # change z
            time.sleep(0.5)
        else:
            print('z change is too big!!!')


def define_wells_manually(meas):
    """ Define positions manually by moving the stage and confirming each selection
    """
    count_limit = 20

    print('\nSelect positions manually!')
    print('.. Confirm each new position with Enter.')
    print('.. To finish selection, press Enter without moving.')

    well_positions = []
    while len(well_positions) < count_limit:
        text = input(f'>> Move to position {len(well_positions) + 1}:')
        if text == '':
            x_current = meas.parameters('ExpControl/scan/range/coarse_x/g_off')
            y_current = meas.parameters('ExpControl/scan/range/coarse_y/g_off')
            # z_current = meas.parameters('ExpControl/scan/range/coarse_z/g_off')
            pos = [x_current, y_current]
            # well_positions.append(pos)
        print(f'<< .. position {len(well_positions) + 1}: {pos}')
        if len(well_positions) >= 1 and pos == well_positions[-1]:
            break
        well_positions.append(pos)

    print(f'# of selected positions: {len(well_positions)}')

    return well_positions


def define_positions(meas, config_num, xy_well, xy_overlap, xy_no_images, shuffle):
    """
    Set up the list of positions, following experiment attributes and microscope settings.
    It creates n x n grid of images

    :param meas: object of imspector measurement
    :param config_num: consecutive number of selected configuration
    :param xy_overlap: amount of frame that will be overlap (if negative, there is no overlap)
    :param xy_no_images: n
    :return: array of [x, y] positions
    """
    conf = meas.configuration(config_num)
    frame_size_x = conf.parameters('ExpControl/scan/range/x/len')
    frame_size_y = conf.parameters('ExpControl/scan/range/y/len')
    # x_start = meas.parameters('ExpControl/scan/range/coarse_x/g_off')
    # y_start = meas.parameters('ExpControl/scan/range/coarse_y/g_off')
    x_start = xy_well[0]
    y_start = xy_well[1]

    # Determine stage xy travel parameters: start, step, range (assuming square pixels and images)
    x_step = (1 - xy_overlap) * frame_size_x
    y_step = (1 - xy_overlap) * frame_size_y
    positions = []

    y_list = list(np.arange(0, xy_no_images, 1) * y_step + y_start)

    # Generate positions for each row
    for y in y_list:
        # Generate the x coordinates
        x_list = list(np.arange(0, - xy_no_images, -1) * x_step + x_start)

        # Run every other row backwards to minimase stage movement
        if y_list.index(y) % 2 == 1:
            x_list.reverse()

        # Populate the final list
        for x in x_list:
            # positions.append([x, y, z])
            positions.append([x, y])

    if shuffle:
        np.random.shuffle(positions)

    return positions


def plot_img_with_rectangles(img_data, coords, save_dir, save_name):
    """ Dummy image analysis routine to identify regions of interest (ROIs) for further measurements.

    Args:
        stack (numpy array): image dataset as provided by im.meas.config.stack
        coords (list of int): list of coordinates in pixels, each as [x,y,width,height]

    Returns:
    """
    my_dpi = 100.
    img = Image.fromarray(img_data)
    fig = plt.figure(figsize=(float(img.size[0]) / my_dpi, float(img.size[1]) / my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.imshow(img)

    for i, [x, y, w, h] in enumerate(coords):
        rect = patches.Rectangle((x, y), w, h, linewidth=.5, edgecolor='r', facecolor='none')
        ax.add_patch(rect)
        ax.text(x + w / 4, y + h / 4, f'{i}', color='r', ha='center', va='center', fontsize=10)

    # plt.show()
    fig.savefig(f'{save_dir}\\{save_name}_objects.png', dpi=800)
    plt.close(fig)

    return None


def move_and_run(meas, config_num, new_position, new_z_pos=None):
    conf = meas.configuration(config_num)
    meas.activate(conf)
    # val = input("Check pos")
    if new_z_pos == None:
        conf.set_parameters('ExpControl/scan/range/x/off', new_position[1])
        time.sleep(0.2)
        conf.set_parameters('ExpControl/scan/range/y/off', new_position[0])
        time.sleep(0.2)
    else:
        conf.set_parameters('ExpControl/scan/range/x/off', new_position[1])
        time.sleep(0.2)
        conf.set_parameters('ExpControl/scan/range/y/off', new_position[0])
        time.sleep(0.2)
        conf.set_parameters('ExpControl/scan/range/z/off', new_z_pos)
        time.sleep(0.2)
    im.run(meas)
    return None


def analyse_xy(im, meas, config_num, nucl_stack_numb, object_stack_numb, n_threshold, ob_threshold, img_part, cellpose_seg=False):
    conf = meas.configuration(config_num)
    meas.activate(conf)
    im.run(meas)
    object_stack = conf.stack(object_stack_numb)
    nucl_stack = conf.stack(nucl_stack_numb)

    data_object = object_stack.data()[0][0]
    data_nucleus = nucl_stack.data()[0][0]

    if cellpose_seg:
        output, nucl_mask = cellpose_segmentation.detect_instances(data=np.array([data_nucleus]), model_type='cyto2', diameter=90,
                                                                   output='segmentation', img_type='img_arr')
        nucl_mask[nucl_mask > 0] = 1
        nucl_mask = np.array(nucl_mask, dtype=bool)
    else:
        nucl_mask = get_nucl_mask(data_nucleus, n_threshold)
    # obj_mask = get_objects(data_object, nucl_mask, ob_threshold)
    # obj_centers = get_object_positions(obj_mask, img_part)
    obj_centers = get_nucl_cent(nucl_mask, img_part)
    return obj_centers, data_nucleus


def get_max_signal(data):
    data_sum = np.sum(data, axis=1)
    max_pos = np.argmax(data_sum)
    return max_pos


def get_nucl_mask(nucl_data, threshold):
    ret, bw_img = cv2.threshold(nucl_data, threshold, 1, cv2.IMREAD_GRAYSCALE)
    gauss_blur = cv2.GaussianBlur(bw_img, (11, 11), 0)
    closing = skimage.morphology.area_closing(gauss_blur, area_threshold=100)
    final = skimage.morphology.binary_dilation(closing)
    return final


def get_objects(objec_data, mask, threshold):
    new_stack = objec_data * mask
    ret2, bw_img2 = cv2.threshold(new_stack, threshold, 1, cv2.IMREAD_GRAYSCALE)
    gauss_blur = cv2.GaussianBlur(bw_img2, (3, 3), 0)
    labels = skimage.measure.label(gauss_blur)
    return labels


def get_object_positions(masks, img_part):
    """
    Calculates centers of objects and their maximum and minimum y and x coordinates from a given mask
    :param masks: mask of objects where each object is represented with a unique number
    :return: arrays of max, min and center coordinates of cells [y, x]
    """
    uniq = np.unique(masks)
    uniq = uniq[uniq > 0]
    dim = np.array(np.shape(masks))
    min_img = dim * img_part
    max_img = (1-img_part) * dim
    all_centers = []
    all_min = []
    all_max = []

    for elem in uniq:
        select_arr = np.copy(masks)
        select_arr[select_arr != elem] = 0
        center = ndimage.measurements.center_of_mass(select_arr)
        center_arr = np.array(center)
        if all(center_arr > min_img) and all(center_arr < max_img):
            all_centers.append(center_arr)
            ii = np.where(masks == elem)
            min_dim = [np.min(ii[0]), np.min(ii[1])]
            max_dim = [np.max(ii[0]), np.max(ii[1])]
            all_min.append(min_dim)
            all_max.append(max_dim)

    return np.array(all_centers)


def get_new_pos(meas, config_num, areas):
    """
    Calculates new offsets of the microscope stage
    :param meas: object of imspector measurement
    :param config_num: consecutive number of selected configuration
    :param sel_area: array of regions of interest with coordinates relative to an image
    :return: array of positions
    """
    conf = meas.configuration(config_num)
    meas.activate(conf)
    sel_area = np.copy(areas)
    psx = conf.parameters('ExpControl/scan/range/x/psz')
    psy = conf.parameters('ExpControl/scan/range/y/psz')
    lenx = conf.parameters('ExpControl/scan/range/x/res') * psx
    leny = conf.parameters('ExpControl/scan/range/y/res') * psy
    if len(sel_area) == 0:
        offsets = []
        return offsets
    else:
        sel_area[:, 0] = sel_area[:, 0] * psy
        sel_area[:, 1] = sel_area[:, 1] * psx
        sel_area[:, 0] = sel_area[:, 0] - leny / 2
        sel_area[:, 1] = sel_area[:, 1] - lenx / 2
        return sel_area


if __name__ == '__main__':

    ##############################
    #      SET PARAMETERS        #
    ##############################

    # set save parameters
    save_dir = "D:\\2024_07_30 - Bostjan, Anja - KI samples RNA\\HCR FISH_New lasers"
    exp_name = 'ITS2_B3 vs 45s_B4'

    # set imspector parameters, counting starts with 0!!!
    meas_n = 0
    af_config = 3                                       # xz measurement for AF
    overview_config = 1                                 # where position of objects is determined
    nucl_stack = 0
    object_stack = 2
    xy_STED_config = 2
    xz_STED_config = 0                                  # where xy is recorded (nucleus, membrane, nanoparticles an bf)
    small_con_config = 4

    # set position parameters
    max_z = 0.00325                                  # in meters
    xy_no_images = 10                                # n for n x n grid of positions
    xy_overlap = -0.6                                 # how much will positions overlap, can be negative, can be decimal
    shuffle_positions = True                       # True - shuffle positions, False - microscope will move in a snake shape

    # analysis parameters
    part_of_image = 0.05                            # how much of image is not considered for analysis (on all sides)
    analyse_with_cellpose = False                             # if True cellpose will be used for nucleus segmentation
    random_selection = False                                 # if True, objects will be selected at random (each time xy and AF scan)
    max_num_per_frame = 10
    nucl_threshold = 10
    object_threshold = 25

    max_count = 50
    ##############################
    #         ALL DONE           #
    ##############################


    #new folde
    notused = 'notused'
    notused_dir = os.path.join(save_dir, notused)
    try:
        os.mkdir(notused_dir)
    except OSError as error:
        print(error)

    # connect to imspector
    im = sp.get_application()
    mnames = im.measurement_names()
    meas = im.measurement(mnames[meas_n])


    # get starting z position
    zrow_start, z_start = get_initial_z(meas, af_config)
    # positions = define_positions(meas, overview_config, xy_overlap, xy_no_images, shuffle_positions)
    start_time = time.time()

    # get all well positions
    well_pos = define_wells_manually(meas)

    for l, well in enumerate(well_pos):
        count = 0
        print(f'--- on well {l + 1} out of {len(well_pos)} ---\n')
        positions = define_positions(meas, overview_config, well, xy_overlap, xy_no_images, shuffle_positions)

        for j, pos in enumerate(positions):
            print(f"Number of objects {count}\{max_count}")
            if count>max_count:
                break
            print(f'----- on position {j+1} out of {xy_no_images*xy_no_images} -----\n')
            meas.set_parameters('ExpControl/scan/range/coarse_x/g_off', pos[0])
            time.sleep(0.2)
            meas.set_parameters('ExpControl/scan/range/coarse_y/g_off', pos[1])
            time.sleep(0.2)
            adjust_focus(meas, af_config, z_start, zrow_start, max_z)
            name = exp_name + f'_well{l:03}_pos{j:03}'
            all_objects, obj_image = analyse_xy(im, meas, overview_config, nucl_stack, object_stack, nucl_threshold,
                                                object_threshold, part_of_image, analyse_with_cellpose)
            rect = [[elem[1] - 20, elem[0] - 20, 40, 40] for elem in all_objects]

            if len(all_objects) == 0:
                meas.save_as(f'{notused_dir}\\{name}.msr')
                plot_img_with_rectangles(obj_image, rect, notused_dir, name)
            else:
                plot_img_with_rectangles(obj_image, rect, save_dir, name)
            if random_selection:
                num_of_obj = len(all_objects)
                if num_of_obj == 1: num_of_obj = 1.1
            else:
                num_of_obj = 1.1

            for i in range(round(num_of_obj*0.5)):
                if random_selection:
                    if i > 0:
                        adjust_focus(meas, af_config, z_start, zrow_start, max_z)
                        all_objects, obj_image = analyse_xy(im, meas, overview_config, nucl_stack, object_stack, nucl_threshold, object_threshold,
                                                            part_of_image, analyse_with_cellpose)
                    name_new = name + f'_obj{i:02}'
                    num_of_obj = len(all_objects)
                    rand_select = np.random.randint(0, num_of_obj)
                    selec_obj = all_objects[rand_select]
                    new_positions = get_new_pos(meas, overview_config, np.array([selec_obj]))
                else:
                    new_positions = get_new_pos(meas, overview_config, all_objects)

                for num, offset_pos in enumerate(new_positions):
                    if num >= max_num_per_frame:
                        break
                    if not random_selection:
                        name_new = name + f'_obj{num:02}'
                    count= count+1
                    # run xz STED
                    move_and_run(meas, xz_STED_config, offset_pos)
                    conf_xz = meas.configuration(xz_STED_config)
                    # get new z
                    stack_xz = conf_xz.stack(0)
                    z_offset = get_max_signal(stack_xz.data()[0][0])
                    new_z = z_offset - conf_xz.parameters('ExpControl/scan/range/z/res')/2
                    new_z_off = new_z * conf_xz.parameters('ExpControl/scan/range/z/psz')
                    # run xy STED
                    move_and_run(meas, xy_STED_config, offset_pos, new_z_off)

                    # run small conf
                    move_and_run(meas, small_con_config, offset_pos, new_z_off)
                    # get image with marked object
                    if not random_selection:
                        selec_obj = all_objects[num]
                    # rect = [[selec_obj[1] - 7, selec_obj[0] - 7, 14, 14]]
                    # plot_img_with_rectangles(obj_image, rect, save_dir, name_new)

                    meas.save_as(f'{save_dir}\\{name_new}.msr')
    
    print('---- all done ----\n')
