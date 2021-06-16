import csv

file = open("WA_Fn-UseC_-HR-Employee-Attrition.csv", mode='r', encoding='utf-8-sig')
#encoding removes the \uefeff characters :)

csvreader = csv.reader(file)

rowlist = [row for row in csvreader]
for category in rowlist[0]:
    index = rowlist[0].index(category)
    exec(category + " = [row[index] for row in rowlist[1:]]")
#if your csv file has categories Age and Gender, now you'll have lists Age and Gender with the appropriate values!

file.close()