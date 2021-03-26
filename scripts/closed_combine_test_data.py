import csv
from copy import deepcopy
from math import pi

def rearrange_data_():
    # f1 = open('./data/CLOSED_CALIBRATION/input_data/q_diff_dir1.txt','r',newline='')
    # f2 = open('./data/CLOSED_CALIBRATION/input_data/q_diff_dir2.txt','r',newline='')
    # f3 = open('./data/CLOSED_CALIBRATION/input_data/q_diff_dir.txt','w',newline='')

    # f1 = open('./data/CLOSED_CALIBRATION/input_data/input_data_1.txt','r',newline='')
    # f2 = open('./data/CLOSED_CALIBRATION/input_data/input_data_2.txt','r',newline='')
    # f3 = open('./data/CLOSED_CALIBRATION/input_data/input_data.txt','w',newline='')

    f1 = open('./data/fixed_calibration/input_data/right_q1_-500_000.txt','r',newline='')
    f2 = open('./data/fixed_calibration/input_data/right_q2_-500_000.txt','r',newline='')
    f3 = open('./data/fixed_calibration/input_data/input_data_right_2.txt','w',newline='')


    data = []
    ## read q
    reader = csv.reader(f1)
    for line in reader:
        data.append(line[0].split())
    f1.close()

    reader = csv.reader(f2)
    for line in reader:
        data.append(line[0].split())
    f2.close()

    for i in range(len(data)):
        for q in range(len(data[i])):
            f3.write(data[i][q])
            if q < (len(data[i])-1):
                f3.write('\t')
        if i < (len(data)-1):
            f3.write('\n')
    f3.close()

    print('total data:',len(data))

rearrange_data_()