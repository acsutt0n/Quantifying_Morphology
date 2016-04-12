# example from python4astronomers.github.io/files/asciifiles.html
# reading input from the web

from BeautifulSoup import BeautifulSoup
def html2tsv(html, index=0):
  """Parse the index'th HTML table in ''html''. Return table as a list of
  tab-separated ASCII table lines"""
  soup = BeautifulSoup(html)
  tables = soup.findAll('table')
  table = tables[index]
  out = []
  for row in table.findAll('tr'):
    colvals = [col.text for col in row.findAll('td')]
    out.append('\t'.join(colvals))
  return out

import urllib2
from astropy.io import ascii
html = urllib2.urlopen('http://hea-www.harvard.edu/XJET/').read()
  # get from web page as string
table1 = html2tsv(html, 0) # parse the first table 
table2 = html2tsv(html, 1) # second table
table3 = html2tsv(html, 2) # third table
