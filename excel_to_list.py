import pandas as pd

#gets all the columns/classes from the excel spreadshet
def excel_to_list(filepath):
    data = pd.read_excel(filepath) 
    df = pd.DataFrame(data)
    array = df.to_numpy()
    data_list = array.tolist()
    return data_list
	
#gets only the columns/classes you want from the excel spreadsheet
def excel_to_list(filepath,columns):
    data = pd.read_excel(filepath) 
    df = pd.DataFrame(data,columns= col) #column names must be identical to the name given in the excel spreadsheet
    array = df.to_numpy()
    data_list = array.tolist()
    return data_list