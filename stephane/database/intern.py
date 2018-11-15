

import people

class Intern(people.People,object):
    # Attributes of people : name, affiliation, email adress 
    
    def __init__(self,d):
        super(Intern,self).__init__(d)    
        
    def display(self):
        super(Intern,self).display()
        
    def display_cosupervisor(self,typ='console'):
        s= self.cosupervisor
        if len(s)>0:
            if s.find('Lab')+3-len(s)==0:
                s = 'in '+s
            else:
                s = 'with '+s
            
            if typ=='console':
                return ', '+s
            if typ=='html':
                return ', <i>'+s +'</i>'
        else:
            return ''
        
    def display_fancy(self,typ='html'):
        return self.surname+' '+self.name+', '+self.title +self.display_cosupervisor(typ=typ) +' ('+self.year+', '+ self.affiliation +')'

    def duration(self):
        self.data_end-self.date_start