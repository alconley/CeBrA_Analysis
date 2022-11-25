from codecs import ascii_decode
import pathlib
import shutil
import xml.etree.ElementTree as ET
import os
import os.path
import argparse

'''
This is a branch of the xml_handle code that is modified to sort through gamma-ray spectra using the CeBrA detectors at the SPS

Takes one xml file that has fitted peaks in HDTV, loops through and creates 2 data files, one with uncalibrated information and one with
calibrated information --> the information is position, position err, width, width err, volume, volume err
'''

def file_handler(dir):
    '''
    File handling portion of the code, where, if wanted, the code can put the xml .txt files in a separate analyzed directory after the code has been run
    on them, as well as outputting and uncalibrated and calibrated data to their own directory. NOT NEEDED
    '''
    os.chdir('../')
    analyzed_folder = os.getcwd()+'/analyzed_fits'
    data_folder = os.getcwd()+'/data_files'
    if not (os.path.isdir(analyzed_folder)):
        os.mkdir(analyzed_folder)
    if not (os.path.isdir(data_folder)):
        os.mkdir(data_folder)
    
    os.chdir(dir)
    for file in os.listdir(dir):
        src_dir = dir + '/' + file
        if file.endswith('.py'):
            continue
        elif file.startswith('detector'):
            dest_dir = data_folder + '/' + file
            shutil.move(src_dir, dest_dir)
        else:
            dest_dir = analyzed_folder + '/' + file
            shutil.move(src_dir, dest_dir)
            # field_setting_folder = analyzed_folder + '/' + B_field + 'G_setting'
            # if not (os.path.isdir(field_setting_folder)):
            #     os.mkdir(field_setting_folder)
            
            # sub_analyzed_folder = field_setting_folder + '/' + angle + '_deg'
            # if not (os.path.isdir(sub_analyzed_folder)):
            #     os.mkdir(sub_analyzed_folder)
            # dest_dir = sub_analyzed_folder + '/' + file
            # shutil.move(src_dir, dest_dir)
    print("\n")
    print("Fit files moved to analyzed_fits directory")
    print("Calibrated & uncalibrated files moved to data_files directory")


def extract_data(B_field, angle, directory):
    '''
    Function that works best if you have multiple fit files over an entire spectra, does not loop over the /region node in the xml tree,
    therefore multiple fit files are needed. This is the original version of the code.
    '''
    #finds the uncalibrated fit & corresponding data (error in peak pos, volume, width)
    fit_list = []
    fit_err_list = []
    width_list = []
    width_err_list = []
    volume_list = []
    volume_err_list = []
    fit_counter = 0
    for file in os.listdir(directory):
        if file.endswith('.py'):
            continue
        elif file.startswith('data'):
            continue
        else:
            fit_counter += 1
            print('Starting fit:', fit_counter, 'now!' )
            mytree = ET.parse(file)
            myroot = mytree.getroot()
            for i in myroot[0]:
                if i.tag == 'peak':
                    for child in i.iter():
                        if child.tag == 'cal':
                            break
                        if child.tag == 'pos':
                            for newchild in child.iter():
                                if newchild.tag == 'value':
                                    fit_value = newchild.text
                                    fit_list.append(float(fit_value))
                                elif newchild.tag == 'error':
                                    fit_err = newchild.text
                                    fit_err_list.append(float(fit_err))
                        elif child.tag == 'vol':
                            for newchild in child.iter():
                                if newchild.tag == 'value':
                                    vol_value = newchild.text
                                    volume_list.append(float(vol_value))
                                elif newchild.tag == 'error':
                                    vol_err = newchild.text
                                    volume_err_list.append(float(vol_err))
                        elif child.tag == 'width':
                            for newchild in child.iter():
                                if newchild.tag == 'value':
                                    width_value = newchild.text
                                    width_list.append(float(width_value))
                                elif newchild.tag == 'error':
                                    width_err = newchild.text
                                    width_err_list.append(float(width_err))
    # print('Peaks', fit_list)
    # print('Peaks Err', fit_err_list)
    # print('Widths', width_list)
    # print('Width Err', width_err_list)
    # print('Volume',volume_list)
    # print('Volume err',volume_err_list)
    
    data = open("data_" + B_field + "G_" + angle + "deg.txt", 'w')
    list = []
    for val in zip(fit_list, fit_err_list, width_list, width_err_list, volume_list, volume_err_list):  #interleaves lists together
        list.append(val)
    print('Data extracted, writing to data_' + B_field + 'G_' + angle + 'deg.txt file!')
    print('Reminder, format for outputted data is: \n', 'Channel  (Err)  Width  (Err)  Volume  (Err)')
    # print(list) --> Unsorted list

    sorted_list = sorted(list, reverse=True)
    # print("New list \n", sorted_list)  ---> Sorted list

    for t in sorted_list:
        line = ' '.join(str(x) for x in t)
        data.write(line + '\n')
    data.close() 
    return data

# just a test function to see if able to properly get positions of peaks
def get_positions(directory):

    cal_data = []
    uncal_data = []

    for file in os.listdir(directory):
       print(file)
       if file.endswith('.py'):
           continue
       else:
           mytree = ET.parse(file)
           myroot = mytree.getroot()
           print(myroot)

           #finds the calibrated & uncalibrated positions of the fits
           for i in myroot[0]:
                if i.tag == 'peakMarker':
                    for child in i.iter():
                        if child.tag == 'cal':
                            cal = child.text
                        elif child.tag == 'uncal':
                            uncal = child.text
                    cal_data.append(cal)
                    uncal_data.append(uncal)

    print(cal_data)
    print(uncal_data)

# if you want user input, but we run with automation so not used as of now
def get_input():
    '''
    Function used to help name files based on what angle/magnetic field setting that is being analyzed.
    '''
    directory = os.getcwd()
    print("Enter B-Field and Angle settings")
    while True:
        try:
            B_field = int(input("B-Field (G): "))
        except ValueError:
            print("Invalid Input, try again")
        else:
            break
    while True:
        try:
            angle = int(input("Angle (degrees): "))
        except ValueError:
            print("Invalid Input, try again")
        else:
            break
    
    return directory, str(B_field), str(angle)

def gamma_Data(file):
    '''
    Function that is used to loop over many regions in one big xml fit file, written for the purpose of analyzing our gamme spectra
    for our CeBrA setup.
    '''
    namesplit = file.split('_')
    detectorID = namesplit[3]
    print(detectorID)

    mytree = ET.parse(file)
    myroot = mytree.getroot()
 
    uncal_fit_list = []
    uncal_fit_err_list = []
    uncal_width_list = []
    uncal_width_err_list = []
    uncal_volume_list = []
    uncal_volume_err_list = []

    cal_fit_list = []
    cal_fit_err_list = []
    cal_width_list = []
    cal_width_err_list = []
    cal_volume_list = []
    cal_volume_err_list = []

    for fit in myroot:
        for i in fit:
            if i.tag == 'peak':
                for child in i.iter():
                    if child.tag == 'uncal':
                        for j in child.iter():
                            if j.tag == 'pos':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        fit_value = newchild.text
                                        uncal_fit_list.append(float(fit_value))
                                    elif newchild.tag == 'error':
                                        fit_err = newchild.text
                                        uncal_fit_err_list.append(float(fit_err))
                            elif j.tag == 'vol':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        vol_value = newchild.text
                                        uncal_volume_list.append(float(vol_value))
                                    elif newchild.tag == 'error':
                                        vol_err = newchild.text
                                        uncal_volume_err_list.append(float(vol_err))
                            elif j.tag == 'width':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        width_value = newchild.text
                                        uncal_width_list.append(float(width_value))
                                    elif newchild.tag == 'error':
                                        width_err = newchild.text
                                        uncal_width_err_list.append(float(width_err))

                    #gets the calibrated data information                
                    if child.tag == 'cal':
                        for j in child.iter():
                            if j.tag == 'pos':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        fit_value = newchild.text
                                        cal_fit_list.append(float(fit_value))
                                    elif newchild.tag == 'error':
                                        fit_err = newchild.text
                                        cal_fit_err_list.append(float(fit_err))
                            elif j.tag == 'vol':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        vol_value = newchild.text
                                        cal_volume_list.append(float(vol_value))
                                    elif newchild.tag == 'error':
                                        vol_err = newchild.text
                                        cal_volume_err_list.append(float(vol_err))
                            elif j.tag == 'width':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        width_value = newchild.text
                                        cal_width_list.append(float(width_value))
                                    elif newchild.tag == 'error':
                                        width_err = newchild.text
                                        cal_width_err_list.append(float(width_err))
    
    # uncalibrated data handling
    uncal_list = []
    for val in zip(uncal_fit_list, uncal_fit_err_list, uncal_width_list, uncal_width_err_list, uncal_volume_list, uncal_volume_err_list):  #interleaves lists together
        uncal_list.append(val)
    sorted_uncal_list = sorted(uncal_list, reverse=True)
    #print(sorted_uncal_list)

    uncal_data = open("detector_" + detectorID + "_uncalibrated_data.txt", 'w')
    for t in sorted_uncal_list:
        line = ' '.join(str(x) for x in t)
        uncal_data.write(line + '\n')
    uncal_data.close() 


    # calibrated data handling

    cal_list = []
    for val in zip(cal_fit_list, cal_fit_err_list, cal_width_list, cal_width_err_list, cal_volume_list, cal_volume_err_list):  #interleaves lists together
        cal_list.append(val)
    sorted_cal_list = sorted(cal_list, reverse=True)

    cal_data = open("detector_" + detectorID + "_calibrated_data.txt", 'w')
    for t in sorted_cal_list:
        line = ' '.join(str(x) for x in t)
        cal_data.write(line + '\n')
    cal_data.close() 

    return uncal_data, cal_data

def general_xml(file, fname):
    mytree = ET.parse(file)
    myroot = mytree.getroot()
 
    uncal_fit_list = []
    uncal_fit_err_list = []
    uncal_width_list = []
    uncal_width_err_list = []
    uncal_volume_list = []
    uncal_volume_err_list = []

    cal_fit_list = []
    cal_fit_err_list = []
    cal_width_list = []
    cal_width_err_list = []
    cal_volume_list = []
    cal_volume_err_list = []

    for fit in myroot:
        for i in fit:
            if i.tag == 'peak':
                for child in i.iter():
                    if child.tag == 'uncal':
                        for j in child.iter():
                            if j.tag == 'pos':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        fit_value = newchild.text
                                        uncal_fit_list.append(float(fit_value))
                                    elif newchild.tag == 'error':
                                        fit_err = newchild.text
                                        uncal_fit_err_list.append(float(fit_err))
                            elif j.tag == 'vol':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        vol_value = newchild.text
                                        uncal_volume_list.append(float(vol_value))
                                    elif newchild.tag == 'error':
                                        vol_err = newchild.text
                                        uncal_volume_err_list.append(float(vol_err))
                            elif j.tag == 'width':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        width_value = newchild.text
                                        uncal_width_list.append(float(width_value))
                                    elif newchild.tag == 'error':
                                        width_err = newchild.text
                                        uncal_width_err_list.append(float(width_err))

                    #gets the calibrated data information                
                    if child.tag == 'cal':
                        for j in child.iter():
                            if j.tag == 'pos':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        fit_value = newchild.text
                                        cal_fit_list.append(float(fit_value))
                                    elif newchild.tag == 'error':
                                        fit_err = newchild.text
                                        cal_fit_err_list.append(float(fit_err))
                            elif j.tag == 'vol':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        vol_value = newchild.text
                                        cal_volume_list.append(float(vol_value))
                                    elif newchild.tag == 'error':
                                        vol_err = newchild.text
                                        cal_volume_err_list.append(float(vol_err))
                            elif j.tag == 'width':
                                for newchild in j.iter():
                                    if newchild.tag == 'value':
                                        width_value = newchild.text
                                        cal_width_list.append(float(width_value))
                                    elif newchild.tag == 'error':
                                        width_err = newchild.text
                                        cal_width_err_list.append(float(width_err))
    
    # uncalibrated data handling
    uncal_list = []
    for val in zip(uncal_fit_list, uncal_fit_err_list, uncal_width_list, uncal_width_err_list, uncal_volume_list, uncal_volume_err_list):  #interleaves lists together
        uncal_list.append(val)
    sorted_uncal_list = sorted(uncal_list, reverse=True)
    #print(sorted_uncal_list)

    uncal_data = open(fname + "_uncalibrated_data.txt", 'w')
    for t in sorted_uncal_list:
        line = ' '.join(str(x) for x in t)
        uncal_data.write(line + '\n')
    uncal_data.close() 


    # calibrated data handling

    cal_list = []
    for val in zip(cal_fit_list, cal_fit_err_list, cal_width_list, cal_width_err_list, cal_volume_list, cal_volume_err_list):  #interleaves lists together
        cal_list.append(val)
    sorted_cal_list = sorted(cal_list, reverse=True)

    cal_data = open(fname + "_calibrated_data.txt", 'w')
    for t in sorted_cal_list:
        line = ' '.join(str(x) for x in t)
        cal_data.write(line + '\n')
    cal_data.close() 

    print("made it here")
    return uncal_data, cal_data

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("xml", type=pathlib.Path, help="Xml file for data to be extracted from")
    args = parser.parse_args()
    return args

def main():
    '''
    Current version will output an uncalibrated & calibrated data set of the given file as a positional argument to the script
    i.e. python3 xml_handle_cebra.py filename.txt is how to run and will return 2 files:
    filename_uncalibrated_data.txt
    filename_calibrated_data.txt
    '''
    dir = os.getcwd()
    args = parseArgs()
    file = dir + '/' + str(args.xml)
    path_name, ext = os.path.splitext(file)
    fname = path_name.split('/')[-1]
    general_xml(file, str(fname))


    '''
    stuff below here will be re-implented to make a more complete, general program
    '''
    # os.chdir(dir + '/fits')
    # fits_dir = os.getcwd()
    # for file in os.listdir():
    #     uncalibrated_data, calibrated_data = gamma_Data(file)
    #     print(uncalibrated_data)
    #     print(calibrated_data)
    # file_handler(fits_dir)  # need to rewrite the file handler to put the files away correctly

    
    #file_handler(dir, data_file, angle, B_field)

if __name__ == '__main__':
    main()

#dir, B_field, angle = get_input()
#get_positions(dir)
#data_file = extract_data(B_field, angle, dir)