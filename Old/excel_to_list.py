import pandas as pd

def excel_to_list_1(filepath):
    """
    Acquires all of the columns/classes from the excel spreadsheet found at 'filepath'.
    'filepath' must be in the format r"(filepath)", ex. filepath = r"C:\Users\Tim\Downloads\HR_Dataset.xlsx"
    """
    data = pd.read_excel(filepath) 
    df = pd.DataFrame(data)
    array = df.to_numpy()
    data_list = array.tolist()
    return data_list

def excel_to_list_2(filepath,columns):
    """
    Acquires only the requested columns/classes from the excel spreadsheet found at 'filepath'.
    'filepath' must be in the format described above. 'columns' is the list of requested columns.
    """
    data = pd.read_excel(filepath) 
    df = pd.DataFrame(data,columns= columns) #column names must be identical to the name given in the excel spreadsheet
    array = df.to_numpy()
    data_list = array.tolist()
    return data_list
