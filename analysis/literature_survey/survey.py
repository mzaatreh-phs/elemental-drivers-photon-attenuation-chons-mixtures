#!/usr/bin/env python3
"""CrossRef survey of low-Z / organic / biological / tissue-equivalent photon-attenuation studies (2015-2026). Anonymous requests, no personal data."""
import json, re, time, urllib.request, urllib.parse, collections
Q=['gamma-ray shielding low-Z materials','photon attenuation biological materials effective atomic number','tissue equivalent materials photon attenuation mass attenuation coefficient',
   'amino acids photon attenuation mass attenuation coefficient effective atomic number','carbohydrates fatty acids proteins gamma ray attenuation parameters','bio-based materials radiation shielding',
   'polymers photon attenuation parameters Geant4 low atomic number','hydrogen-rich materials gamma attenuation effective atomic number','organic compounds gamma-ray shielding mass attenuation coefficient',
   'lead-free polymer composite gamma shielding review','radiation shielding materials review polymer composites','elemental composition photon interaction effective atomic number mixture rule',
   'biomolecules photon interaction parameters buildup factor','lightweight gamma shielding natural materials wood plant food attenuation','human tissues Geant4 gamma-ray transmission buildup factor',
   'sulfur containing compounds gamma attenuation','low atomic number materials photon shielding photoelectric Compton composition','phantom materials photon attenuation water PMMA effective atomic number']
def get(u):
    r=urllib.request.Request(u,headers={'User-Agent':'literature-survey/1.0'}); return json.load(urllib.request.urlopen(r,timeout=30))
recs={}
for q in Q:
    try:
        u='https://api.crossref.org/works?'+urllib.parse.urlencode({'query.bibliographic':q,'rows':50,'filter':'from-pub-date:2015-01-01,type:journal-article','select':'DOI,title,issued,container-title,abstract,author'})
        for it in get(u)['message']['items']:
            t=(it.get('title') or [''])[0]
            if t and it['DOI'] not in recs: recs[it['DOI']]={'title':t,'year':(it.get('issued',{}).get('date-parts') or [[None]])[0][0],'journal':(it.get('container-title') or [''])[0],'abstract':re.sub('<[^>]+>',' ',it.get('abstract','') or ''),'query':q,'first':(it.get('author') or [{}])[0].get('family','')}
    except Exception as e: print('ERR',q,e)
    time.sleep(0.4)
json.dump(recs,open('crossref_records.json','w'),indent=1); print(len(recs),'unique records; with abstract:',sum(1 for r in recs.values() if len(r['abstract'])>80))
