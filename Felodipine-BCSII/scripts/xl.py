import zipfile,re,sys,xml.etree.ElementTree as ET
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def read(path):
    z=zipfile.ZipFile(path)
    ss=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns):
            ss.append(''.join(t.text or '' for t in si.iter('{%s}t'%ns['m'])))
    wb=ET.fromstring(z.read('xl/workbook.xml'))
    rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    rmap={r.get('Id'):r.get('Target') for r in rels}
    out={}
    for sh in wb.find('m:sheets',ns):
        rid=sh.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        t=rmap[rid].lstrip('/'); t=t if t.startswith('xl/') else 'xl/'+t
        rows=[]
        for r in ET.fromstring(z.read(t)).iter('{%s}row'%ns['m']):
            row=[]
            for c in r.findall('m:c',ns):
                v=c.find('m:v',ns); typ=c.get('t')
                if typ=='s': row.append(ss[int(v.text)])
                elif typ=='inlineStr': row.append(''.join(x.text or '' for x in c.iter('{%s}t'%ns['m'])))
                elif v is not None:
                    try: row.append(float(v.text))
                    except: row.append(v.text)
                else: row.append(None)
            rows.append(row)
        out[sh.get('name')]=rows
    return out
