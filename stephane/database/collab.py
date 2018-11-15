

import people

class Collab(people.People,object):
    # Attributes of people : name, affiliation, email adress 
    
    def __init__(self,d):
        super(Collab,self).__init__(d)    
        
    def display(self):
        super(Collab,self).display()
        
    def display_fancy(self):
        return self.surname+' '+self.middle+''+self.name+' ('+self.affiliation + ', '+self.department+')'
    
    def duration(self):
        self.data_end-self.date_start