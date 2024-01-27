import requests
import pandas as pd
from multiprocessing import Pool
from bs4 import BeautifulSoup

base_URL = 'https://clustercad.jbei.org/pks/'

def get_data(bgc_values):
    thread_list = []
    bgc_id, is_reviewed = bgc_values
    s = requests.session()
    s.headers = {'X-Requested-With': 'XMLHttpRequest'} #required header to not get a 404

    bgc_url = base_URL + str(bgc_id) + '/'
    bgc_page = BeautifulSoup(s.get(bgc_url).content, "html.parser")

    genbank_id = bgc_page.find_all('a')[7].text.strip()

    final_smiles = bgc_page.find_all('img')[1]['data-smiles']

    cluster_notes = bgc_page.find_all('dd')[4].text.strip()

    for gene in bgc_page.find_all(class_='list-group'):
        locus_tag = gene['id']

        module_container = gene.find(class_='row')
        for module in module_container.findChildren('div', recursive=False):
            module_id = module.text.strip().split('\n')[0].split(' ')[1]
            module_smiles = module.find('img')['data-smiles']
            for domain in module.find_all('button'):
                name = domain.text.strip()
                name = name if name == 'ACP' else 'PKS_' + name
                
                domain_id = domain['data-domainid']
                domain_URL = base_URL + 'domainLookup?domainid=' + domain_id
                domain_info = s.get(domain_URL).json()
                
                translation = domain_info['AAsequence']
                start = domain_info['start']
                end = domain_info['stop']
                seqlen = len(translation)

                is_active = 'active'
                at_substrate = ''
                kr_type = ''

                if domain.has_attr('title'):
                    domain_description = domain['title']
                    is_active = 'inactive' if 'inactive' in domain_description else is_active
                    if name == 'PKS_AT':
                        at_substrate = domain_description.split(',')[0].split()[1]
                    if name == 'PKS_KR':
                        kr_type = domain_description.split(',')[0].split()[1]

                thread_list.append([name, translation, seqlen, start, end, module_id, bgc_id, genbank_id, is_reviewed, final_smiles, module_smiles, is_active, at_substrate, kr_type, cluster_notes])
    return thread_list

if __name__ == '__main__':
    main_s = requests.session()
    all_URL = base_URL + "all"
    page = BeautifulSoup(main_s.get(all_URL).content, "html.parser")

    table = page.find(id='clusterTable')
    tableBody = table.find_all('tbody')[0]

    bgc_values = []
    for row in tableBody.find_all('tr'):
        #constants values
        row_values = row.find_all('td')
        is_reviewed = True if row_values[-1].text.strip() == 'Y' else False
        bgc_id = row_values[0].text.strip()
        
        bgc_values.append([bgc_id, is_reviewed])  

    bgc_list = []
    with Pool() as p:
        thread_list = p.map(get_data, bgc_values)
        bgc_list = [x for xs in thread_list for x in xs]

    df = pd.DataFrame(bgc_list, columns=['name', 'translation', 'seqlen', 'start', 'end', 'module_id', 'bgc_id', 'genbank_id', 'is_reviewed', 'final_smiles', 'module_smiles', 'is_active', 'at_substrate', 'kr_type', 'cluster_notes'])# (predicted_final_smiles)
    df.to_csv('output.csv', index=False)