import pandas as pd
import numpy as np
import mechanicalsoup
from bs4 import BeautifulSoup as bs

file_dir = r'C:/Users/MSmith/Documents/R/data/'
plants_dir = r'C:/Users/MSmith/Documents/R/data/GBIF_allNG/png_plants.csv'
taxa_dir = r'C:/Users/MSmith/Documents/R/data/GBIF_allNG/PNGgeneracsv.csv'
taxa = pd.read_csv(taxa_dir,header=None)
taxa.columns=['taxon']

for i in taxa['taxon']:
    print(i)
    try:
        myurl = "http://www.pngplants.org/cgi-bin/RBG?form=PNG%2Fpng&laeno=&section=all&fam=&genus={0}&species=&otherran=&otherval=&collname_tab=&collno=&datecol_1=&datecol_2=&datecol_3=&date2_1=&date2_2=&date2_3=&textdate%3A1=&country=&state=&subdiv=&loadloc=&sitedesc=&plantdes=&originat=&recordba=&donat2_tab=".format(i)
        browser = mechanicalsoup.StatefulBrowser()
        browser.open(myurl)
        browser.select_form('form[name="summform"]')
        check = browser.page.find_all('input')
        rangeval = len(check)-2
        names = []
        for j in range(3, rangeval):
            #print(check[i]['value'])
            names.append(check[j]['value'])
        #print(names)
        browser["alpkey"] = names
        #browser.launch_browser()

        response = browser.submit_selected()
        #print(response.url)
        #print(response.text)
        soup = response.soup
        table = soup.find("table")
        #print(f"[+] Found a total of {len(tables)} tables.")
        speciesname = []
        family = []
        date = []
        location = []
        coords = []
        for trs in table.find_all('tr'):
            tds = trs.find_all("td")
            for tds in trs.find_all("td"):
                bseries = tds.find_all("b")
                bsindiv1 = bseries[0] # species name and family
                bstext1 = bsindiv1.text.strip()
                bstext1a = bstext1.split()
                splitlength = len(bstext1a)
                if (splitlength > 2):
                    bstext1b = bstext1a[0] + ' ' + bstext1a[1] # genus and species
                    bstext1c = bstext1a[splitlength-1]
                else:
                    bstext1b = bstext1a[0] + ' no_species' # if genus only
                    bstext1c = bstext1a[1]
                speciesname.append(bstext1b)
                family.append(bstext1c)
                btextfind1 = tds.find('i')
                try:
                    bstext3 = btextfind1.next_sibling.strip()
                except:
                    bstext3 = ''
                date.append(bstext3)
                btextfind2 = tds.find('b', text='Locality: ')
                try:
                    bstext4 = btextfind2.next_sibling.strip()
                except:
                    bstext4 = ''
                location.append(bstext4)
                btextfind3 = tds.find('b', text=':')
                try:
                    bstext5 = btextfind3.next_sibling.strip()
                except:
                    bstext5 = ''
                coords.append(bstext5)
        df = pd.DataFrame(list(zip(speciesname, family, date, location, coords)),
                  columns=['Species', 'Family', 'Location', 'Date', 'Coords'])
        #print(df)
        with open(plants_dir, 'a') as f:
            df.to_csv(f, header=False, index=False, line_terminator='\n')
    except:
        print(i)
        print('- URL not found')