# -*- coding: utf-8 -*-

import people as p
import intern as i
import collab
import unicodecsv as csv

def peoplelist():
    d = {}
    d['name'] = 'Matthieu Labousse'
    d['aff'] = 'ESPCI'
    d['email'] = 'matthieu.labousse@espci.fr'
    
    ex = p.People(d)
    ex.display()
    
    d = {}
    d['name'] = 'Vincent Bacot'
    d['aff'] = 'University Paris-Diderot'
    d['email'] = 'vincent.baco@gmail.com'
    d['dates'] = ''
    
    ex = i.Intern(d)
    ex.display()
    
    print('Students')
    filename = '/Users/stephane/Documents/Etudiants/Recherche/internships.csv' 
    with open(filename) as csvfile:
        reader = list(csv.DictReader(csvfile,delimiter=';'))
        reader.sort(key=lambda x: x['Year'],reverse=True)
        for row in reader:
            student = i.Intern(row);            
            student.display_fancy()        

    print('')
    print('Collaborators')
    filename = '/Users/stephane/Documents/Etudiants/Recherche/collaborators.csv' 
    with open(filename) as csvfile:
        reader = list(csv.DictReader(csvfile,delimiter=';'))
        reader.sort(key=lambda x: x['Name'])
        for row in reader:
            student = collab.Collab(row);            
            print(student.display_fancy())
            

import bibtexparser as bib
import bibstyle
              
def biblio():
    
    filename = '/Users/stephane/Documents/Articles/Published/Liste_biblio/mesarticles.bib'
           
    with open(filename) as bibtex_file:
        bib_database = bib.load(bibtex_file)
        
        entries = list(bib_database.entries)
        entries.sort(key=lambda x: x['year'],reverse=True)

        for entry in entries:
            print(bibstyle.display(entry,typ='console'))

        print('\n \n \n')
        
        # list option
        print('<ol reversed>')
        for entry in entries:
            print(bibstyle.display(entry,typ='html'))
        print('</ol>')
        
import genhtml
from yattag import indent

def peoplelisthtml():
    
    print('Students')
    filename = '/Users/stephane/Documents/Etudiants/Recherche/internships.csv' 
    with open(filename) as csvfile:
        reader = list(csv.DictReader(csvfile,delimiter=';'))
        reader.sort(key=lambda x: x['Year'],reverse=True)
        
        stringlist = []
        for row in reader:
            print(row['Name'])
            student = i.Intern(row);            
            stringlist.append(student.display_fancy()) 
          
        doc, tag, text = genhtml.makelist(stringlist,'Past supervised students',order='u')

        filename='siteweb/students.html'
        f=open(filename,'w')
        f.write(indent(doc.getvalue()).encode('utf-8'))
        f.close()
        
    print('Collaborators')
    filename = '/Users/stephane/Documents/Etudiants/Recherche/collaborators.csv' 
    with open(filename) as csvfile:
        reader = list(csv.DictReader(csvfile,delimiter=';'))
        reader.sort(key=lambda x: x['Name'])
        
        stringlist = []
        for row in reader:
#            print(row['Name'])
            people = collab.Collab(row);            
            stringlist.append(people.display_fancy()) 
          
        doc, tag, text = genhtml.makelist(stringlist,'Collaborators',order='u')

        filename='siteweb/collaborators.html'
        f=open(filename,'w')
        f.write(indent(doc.getvalue()).encode('utf-8'))
        f.close()
        
    print('Biblio')

    filename = '/Users/stephane/Documents/Articles/Published/Liste_biblio/mesarticles.bib'   
    with open(filename) as bibtex_file:
        bib_database = bib.load(bibtex_file)
        
        entries = list(bib_database.entries)
        entries.sort(key=lambda x: x['year'],reverse=True)

        stringlist = []
        for entry in entries:
            stringlist.append(bibstyle.display(entry,typ='html'))

        doc, tag, text = genhtml.makelist(stringlist,'Publications',order='o',opt=' reversed')

        filename='siteweb/publications.html'
        f=open(filename,'w')
        f.write(indent(doc.getvalue()).encode('utf-8'))
        f.close()
                
peoplelisthtml()
#biblio()
#peoplelist()