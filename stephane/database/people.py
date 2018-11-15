



class People(object):
    # Attributes of people : name, affiliation, email adress 
    
    def __init__(self,d):
        setattr(self,'keys',[key.lower() for key in d.keys()]) #key = key.lower()#key.decode('utf-8').lower()            
        for key in d.keys():
            setattr(self,key.lower(),d[key])           

    def display(self):
        for key in self.keys:
            print(key + ' : '+ getattr(self,key))