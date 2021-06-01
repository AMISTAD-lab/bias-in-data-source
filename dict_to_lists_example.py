
dictionary = {'gender':['m','f','nb'], 'gpa':['2.5-3.0','3.0-3.5','3.5-4.0']}

for key in dictionary.keys():
    exec(key+' = '+str(dictionary[key]))
    #if you used the example dictionary provided above,
    #the variable gender == ['m','f','nb']
    #and gpa == ['2.5-3.0','3.0-3.5','3.5-4.0']

#if you wanted a list of the lists in your dictionary, then:
megalist = []
for key in dictionary.keys():
    megalist += [dictionary[key]]
#after running this,
#megalist == [['m','f','nb'], ['2.5-3.0','3.0-3.5','3.5-4.0']]
