


from yattag import Doc

def example():
    doc, tag, text = Doc().tagtext()

    with tag('h1'):
        text('Hello world!')

    print(doc.getvalue())


def makelist(stringlist,title,order='u',opt=''):    
    doc, tag, text = Doc().tagtext()
        
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            doc.asis('<meta charset="UTF-8">')
            with tag('style'):
                text('div.text {line-height: 1.6; font-size: large;}')
        with tag('body'):
            doc.asis('<div class="text">')
            with tag('h2'):
                text(title)
            with tag(order+'l'+opt):
                for elem in stringlist:
                    with tag('li'):
                        print(elem)
                        doc.asis(elem)
            doc.asis('</div>')
            
    return doc, tag, text

#main()