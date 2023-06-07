#!/usr/bin/env python


'''
export the aligned subtomograms selected using GA


~/ln/tomominer/tomominer/average/genetic_algorithm/analysis/export_selected_data.py

'''



import os, json, pickle



def main():
    with open('export_selected_data__op.json') as f:      op = json.load(f)

    with open(op['data file in'], 'rb') as f:       pmpg = pickle.load(f)

    with open(op['data json out'], 'w') as f:       json.dump(pmpg['dj'], f, indent=2)


if __name__ == '__main__':
    main()

