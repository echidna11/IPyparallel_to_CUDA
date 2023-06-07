'''
functions for inserting images into a excel table

~/ln/tomominer/tomominer/doc/table/excel/insert_image.py
'''


import numpy as N



'''
given a list of images and their cell positions, 
get the corresponding row and columes sizes to fit the images
'''

def get_cell_sizes(imgs):

    from scipy.misc import imread
    
    # make a table
    tab = {}
    for img in imgs:
        row = img['row']
        col = img['col']

        if row not in tab:  tab[row] = {}
        assert col not in tab[row]          # make sure each cell has no more than one image

        tab[row][col] = img


    # get original and scaled sizes of images
    for row in tab:
        for col in tab[row]:
            t = tab[row][col]
            im = imread(t['file'])
            t['size'] = N.array(im.shape[:2])

            if 'size_scaled' not in t: 
                t['size_scaled'] = N.copy(t['size'])

                if 'scale' in t:
                    assert 'x_scale' not in t
                    assert 'y_scale' not in t
                    t['x_scale'] = t['scale']
                    t['y_scale'] = t['scale']

                if 'x_scale' in t:  t['size_scaled'][0] *= t['x_scale']
                if 'y_scale' in t:  t['size_scaled'][1] *= t['y_scale']


            else:
                assert 'scale'  not in t
                assert 'x_scale' not in t
                assert 'y_scale' not in t

                # if the final image size has a None entry, then use the other entry to scale the image
                if (t['size_scaled'][0] is None) and (t['size_scaled'][1] is not None):
                    t['size_scaled'][0] = float(t['size'][0] * t['size_scaled'][1]) / t['size'][1]
                elif (t['size_scaled'][0] is not None) and (t['size_scaled'][1] is None):
                    t['size_scaled'][1] = float(t['size'][1] * t['size_scaled'][0]) / t['size'][0]
                else:
                    assert t['size_scaled'][0] is not None
                    assert t['size_scaled'][1] is not None

                t['x_scale'] = float(t['size_scaled'][0]) / t['size'][0]
                t['y_scale'] = float(t['size_scaled'][1]) / t['size'][1]


    cols = []
    for row in tab:     cols.extend(tab[row].keys())

    cols = sorted(list(set(cols)))

    # get row height
    max_heights = {}
    for row in tab:
        for col in tab[row]:
            t = tab[row][col]
            height_t = t['size_scaled'][0]
            if row not in max_heights:      max_heights[row] = height_t
            max_heights[row] = max(max_heights[row], height_t)


    # get col width
    max_widths = {}
    for col in cols:
        for row in tab:
            if col not in tab[row]:     continue

            t = tab[row][col]
            width_t = t['size_scaled'][1]
            if col not in max_widths:         max_widths[col] = width_t
            max_widths[col] = max(max_widths[col], width_t)


    return {'tab':tab, 'max_row_heights':max_heights, 'max_col_widths':max_widths}


'''
adjust row and column sizes
'''
def adjust_row_col_size(worksheet, row_heights, col_widths, row_margin=0.0, col_margin=0.0):
    for row in row_heights:
        worksheet.set_row(row, height=(row_margin + row_heights[row] * 3.0 / 4.0))        # see source code of worksheet._size_row(row)

    max_digit_width = 7  # For Calabri 11.      see source code of worksheet._size_col(col)
    for col in col_widths:
        worksheet.set_column(col, col, width=col_margin + float(col_widths[col]) / max_digit_width)     # see source code of worksheet._size_col(col)


def insert_images(worksheet, imgs, op={'row_margin':0.0, 'col_margin':0.0}):
    c = get_cell_sizes(imgs)

    adjust_row_col_size(worksheet=worksheet, row_heights=c['max_row_heights'], col_widths=c['max_col_widths'], row_margin=op['row_margin'], col_margin=op['col_margin'])

    tab = c['tab']
    for row in tab:
        for col in tab[row]:
            t = tab[row][col]

            x_scale = t['x_scale'] if 'x_scale' in t else 1.0
            y_scale = t['y_scale'] if 'y_scale' in t else 1.0

            x_offset = t['x_offset'] if 'x_offset' in t else 0.0
            y_offset = t['y_offset'] if 'y_offset' in t else 0.0

            worksheet.insert_image(t['row'], t['col'], t['file'], {'x_scale': x_scale, 'y_scale': y_scale, 'x_offset':x_offset, 'y_offset':y_offset})



'''
# test code

import xlsxwriter

workbook = xlsxwriter.Workbook('demo.xlsx')
worksheet = workbook.add_worksheet()

worksheet.write(2, 0, 123)
worksheet.write(3, 0, 123.456)

imgs = []
imgs.append({'file':"python-logo.png", 'row':2, 'col':3})
imgs.append({'file':"python-logo.png", 'row':4, 'col':3})


imgs.append({'file':"python-logo.png", 'row':4, 'col':5, 'x_scale':0.5, 'y_scale':0.8})
imgs.append({'file':"python-logo.png", 'row':5, 'col':5, 'x_scale':0.5, 'y_scale':0.3})

from tomominer.doc.table.excel.insert_image import insert_images
insert_images(worksheet, imgs)

workbook.close()

'''

