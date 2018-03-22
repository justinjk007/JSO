import csv
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

def parse(file_name):
    "Parses the data sets from the csv file we are given to work with"
    file = open(file_name)
    rawFile = csv.reader(file)    # Reading the csv file into a raw form
    rawData = list(rawFile)       # Converting the raw data into list from.
    file.close()
    return rawData

def draw(xCords, yCords, xLabel, yLabel, filename, col):
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title("Pentagon - "+filename)
    plt.plot(xCords, yCords, '-o', color=col, markersize = 1, linewidth = 2)
    # for xy in zip(xCords, yCords):
    #     ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='data')

def toXandY(unorderedData):
    "This method converts seperates x and y co-ordinates for plotting"
    orderedData = []
    orderedData.append([])        # Add a new sublist every time
    orderedData.append([])        # Add a new sublist every time
    listSize = len(unorderedData)
    for x in range(0, listSize):
        orderedData[0].append(unorderedData[x][0])  # Seperates the x-cords
    for y in range(0, listSize):
        orderedData[1].append(unorderedData[y][1])  # Seperates the y-cords
    return orderedData

def main():
    newData = []
    line_x = []
    line_y = []
    file_name = "rawData.csv"
    data = parse(file_name)   # Calling the parse funtion we made
    labels = data.pop(0)            # Necessary evil
    list_size = len(data)
    for i in range(0, list_size):    # Converting the string list to float
        newData.append([])          # Add a new sublsit every time
        for j in range(0, 2):       # Append converted data to the new list
            newData[i].append(float(data[i][j]))
    DataXandY = toXandY(newData)     # DataXandY -> [[Xs][Ys]]
    while len(DataXandY[0]) != 0:
        # Pop two points to form a line
        line_x.append(DataXandY[0].pop(0))
        line_y.append(DataXandY[1].pop(0))
        line_x.append(DataXandY[0].pop(0))
        line_y.append(DataXandY[1].pop(0))
        draw(line_x, line_y, labels[0], labels[1], file_name, "g")
        line_x[:] = []          # Delete the contents
        line_y[:] = []          # Delete the contents

    # plt.xlim(-5, 15)
    # plt.ylim(-5, 15)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

if __name__ == "__main__":
    main()
