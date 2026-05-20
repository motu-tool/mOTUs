import numpy as np
import sqlite3
from typing import List, Tuple, Generator
import zlib
import logging
from rapidfuzz import process, fuzz
import pathlib
from motus import mutils
import tqdm
import urllib

TABLES = {}
TABLES['KEGG'] = ('KEGG_NAME_2_KEGG_ID', 'KEGG_ID_2_GENOME_IDS')
TABLES['TAXONOMY'] = ('TAXMOTU_NAME_2_TAXMOTU_ID', 'TAXMOTU_ID_2_GENOME_IDS')
TABLES['GENOME'] = ('GENOME_NAME_2_GENOME_ID', 'GENOME_ID_2_DATA_IDS')
TABLES['EGGNOG'] = ('EGGNOG_NAME_2_EGGNOG_ID', 'EGGNOG_ID_2_GENOME_IDS')
TABLES['PFAM'] = ('PFAM_NAME_2_PFAM_ID', 'PFAM_ID_2_GENOME_IDS')
GENOME_DB_TYPE = 'GENOME'
EGGNOG_DB_TYPE = 'EGGNOG'
KEGG_DB_TYPE = 'KEGG'
TAXONOMY_DB_TYPE = 'TAXONOMY'
PFAM_DB_TYPE = 'PFAM'


# TODO there will be 0 annotations which cannot be decoded. Have to find an example and fix. Otherwise this break the de-initing routine


class Genome:
    _decoded_kegg_ids = None
    _decoded_taxmotu_ids = None
    _decoded_eggnog_ids = None
    _decoded_pfam_ids = None
    _genome_id = None
    _kegg_names = None
    _taxmotu_names = None
    _eggnog_names = None
    _pfam_names = None
    def __init__(self, genome_id: int):
        self._genome_id = genome_id

    def set_all_ids(self, eggnog_ids, kegg_ids, pfam_ids, taxmotu_ids):
        self._decoded_kegg_ids = kegg_ids
        self._decoded_taxmotu_ids = taxmotu_ids
        self._decoded_eggnog_ids = eggnog_ids
        self._decoded_pfam_ids = pfam_ids
    
    def set_all_names(self, genome_name, eggnog_names, kegg_names, pfam_names, taxmotu_names):
        self._kegg_names = kegg_names
        self._genome_name = genome_name
        self._taxmotu_names = taxmotu_names
        self._eggnog_names = eggnog_names
        self._pfam_names = pfam_names
    
    def get_all_ids(self):
        return self._genome_id, self._decoded_eggnog_ids, self._decoded_kegg_ids, self._decoded_pfam_ids, self._decoded_taxmotu_ids


    def get_genome_name(self):
        if not self._genome_name:
            logging.info('genome name not set.')
            exit(1)           
        return self._genome_name
    def get_eggnog_printable(self):
        if not self._eggnog_names:
            logging.info('eggnog not set.')
            exit(1)
        return ';'.join(self._eggnog_names)


    def get_kegg_printable(self):
        if not self._kegg_names:
            logging.info('kegg not set.')
            exit(1)
        return ';'.join(self._kegg_names)
    
    def get_pfam_printable(self):
        if not self._pfam_names:
            logging.info('pfam not set.')
            exit(1)
        return ';'.join(self._pfam_names)

    def get_tax_printable(self):
        if not self._taxmotu_names:
            logging.info('taxmotu not set.')
            exit(1)
        return ';'.join(self._taxmotu_names[:-1])
    def get_motu_printable(self):
        if not self._taxmotu_names:
            logging.info('taxmotu not set.')
            exit(1)
        return self._taxmotu_names[-1]

class SearchDB:
    _type_2_id_2_name = {}
    _type_2_name_2_id = {}
    _connection: sqlite3.Connection = None

    def __init__(self, database_path: str):
        logging.info('Initializing Search Database')
        self._connection = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True)
        for name_type, (database_table, _) in TABLES.items():
            cursor = self._connection.cursor()
            rows = cursor.execute(f'SELECT name, id FROM {database_table}').fetchall()
            if name_type == 'PFAM':
                rows2 = []
                for row in rows:
                    row2 = (row[0].rsplit('.', 1)[0], row[1])
                    rows2.append(row2)
                rows = rows2
            n2i = dict(rows)
            i2n = {v: k for k, v in n2i.items()}
            self._type_2_id_2_name[name_type] = i2n
            self._type_2_name_2_id[name_type] = n2i
            cursor.close()
        logging.info(f'Initialized Search Database - {sum([len(x) for x in self._type_2_name_2_id.values()]):,} records indexed')

    def close(self):
        self._connection.close()

    def exact_search_keyword(self, keyword:str):
        results = []
        for type, name_2_id in self._type_2_name_2_id.items():
            if keyword in name_2_id:
                results.append((type, keyword, name_2_id[keyword]))
        return results

    def fuzzy_search_keyword(self, keyword: str):
        written = False
        for type, name_2_id in self._type_2_name_2_id.items():
            matches = process.extract(keyword, name_2_id.keys(), scorer=fuzz.WRatio, processor=str.lower,limit=10000000, score_cutoff=80)
            if len(matches) != 0:
                if not written:
                    logging.info(f'Found inexact hits for "{keyword}":')
                    written = True
                logging.info(f'{type}:')
            for (token, score, _) in sorted(matches, key=lambda match: match[1], reverse=True):
                logging.info(f'\tFound "{token}" with score {round(score, 2)}')
        if not written:
            logging.info('No matches found with fuzzy search.')


    _DTYPE_TO_ID = {
        np.uint8:  0,
        np.uint16: 1,
        np.uint32: 2,
        np.uint64: 3,
    }
    _ID_TO_DTYPE = {v: k for k, v in _DTYPE_TO_ID.items()}

    def _delta_decode_np(self, blob: bytes) -> np.ndarray:
        if len(blob) < 2:
            raise ValueError("Corrupt blob: too short")

        complement_flag = blob[0]  # currently always 0
        dtype_id = blob[1]

        if dtype_id not in self._ID_TO_DTYPE:
            raise ValueError(f"Unknown dtype_id {dtype_id}")

        dtype = self._ID_TO_DTYPE[dtype_id]
        payload = blob[2:]

        deltas = np.frombuffer(payload, dtype=dtype)
        ids = np.cumsum(deltas, dtype=np.int64)

        return ids
    def _decode_blob(self, blob: bytes) -> List[int]:
        if not blob:
            raise ValueError("empty blob")
        # flag = blob[0]
        payload = blob
        data = zlib.decompress(payload)    
        return self._delta_decode_np(data).tolist()
  

    def _get_genome_ids_from_db(self, db_name_type, db_name_id):
        cursor = self._connection.cursor()
        t = TABLES[db_name_type][1]
        rows = cursor.execute(f'SELECT id, encoded_data FROM {TABLES[db_name_type][1]} WHERE id = ?', (db_name_id, )).fetchall()
        #print(len(rows), TABLES[db_name_type][1], db_name_id, self._type_2_universe_size[db_name_type])
        genome_ids = []
        for db_name_id, blob in rows:
            _genome_ids = self._decode_blob(blob)
            for _genome_id in _genome_ids:
                genome_ids.append(_genome_id)
        cursor.close()
        return genome_ids

    def get_genome_ids_by_database_hits(self, database_hits: List[Tuple[str, str, int]]):
        # 1. for each hit
        # 2. check if this is a genome hit --> done
        # 3. if not get the correct database and get all genome_ids
        # 4. Return genome_ids.
        # This method will only return genome hits. we might want to process them further so we dont get all data yet
        genome_ids = set()
        for (db_name_type, keyword, db_name_id) in database_hits:
            if db_name_type == GENOME_DB_TYPE:
                genome_ids.add(db_name_id)
                continue
            _genome_ids = self._get_genome_ids_from_db(db_name_type, db_name_id)

            for _genome_id in _genome_ids:
                genome_ids.add(_genome_id)
        return genome_ids

    def get_genomes_by_genome_ids(self, genome_ids: List[int]) -> List[Genome]:
        genome_ids = sorted(list(genome_ids))
        chunked_list = [genome_ids[i:i + 500] for i in range(0, len(genome_ids), 500)]
        genomes = []
        with tqdm.tqdm(total = len(genome_ids), desc="Querying database for annotations", unit = 'genomes') as pbar:
            for chunk in chunked_list:
                cursor = self._connection.cursor()
                placeholders = ','.join(['?'] * len(chunk))
                rows = cursor.execute(f'SELECT id, encoded_eggnog, encoded_kegg, encoded_pfam, encoded_taxmotu FROM {TABLES[GENOME_DB_TYPE][1]} WHERE id IN ({placeholders})', chunk).fetchall()
                cursor.close()

                for genome_id, encoded_eggnog, encoded_kegg, encoded_pfam, encoded_taxmotu in rows:
                    decoded_eggnog = self._decode_blob(encoded_eggnog) 
                    decoded_taxmotu = np.frombuffer(zlib.decompress(encoded_taxmotu), dtype=np.int64).tolist()
                    decoded_kegg = self._decode_blob(encoded_kegg)
                    decoded_pfam = self._decode_blob(encoded_pfam)
                    g = Genome(genome_id)
                    g.set_all_ids(decoded_eggnog, decoded_kegg, decoded_pfam, decoded_taxmotu)
                    genomes.append(g)
                    pbar.update(1)

        return genomes
    
    def get_genome_names(self, genome_ids):
        """Takes a list of genome_ids
        and reports their genome names

        Args:
            genome_ids (List[str]): 

        Returns:
            _type_: _description_
        """        
        genome_names = []
        for genome_id in genome_ids:
            genome_names.append(self._type_2_id_2_name[GENOME_DB_TYPE][genome_id])
        return genome_names


    def stringify_genomes(self, genomes: List[Genome]) -> Generator[Genome, None, None] :
        '''
        Takes a list of genome objects which still have 
        only intified identifiers and will add the 
        stringified versions
        '''
        for genome in genomes:
            genome_id, eggnog_ids, kegg_ids, pfam_ids, taxmotu_ids = genome.get_all_ids()
            genome_name = self._type_2_id_2_name[GENOME_DB_TYPE][genome_id]
            eggnog_names = [self._type_2_id_2_name[EGGNOG_DB_TYPE][x] for x in eggnog_ids if x != 0]
            kegg_names = [self._type_2_id_2_name[KEGG_DB_TYPE][x] for x in kegg_ids if x != 0]
            pfam_names = [self._type_2_id_2_name[PFAM_DB_TYPE][x] for x in pfam_ids if x != 0]
            taxmotu_names = [self._type_2_id_2_name[TAXONOMY_DB_TYPE][x] for x in taxmotu_ids if x != 0]
            genome.set_all_names(genome_name, eggnog_names, kegg_names, pfam_names, taxmotu_names)
            yield genome






def find_genomes(search_tokens: List[str], output_file:pathlib.Path, annotations_to_report: List[str], db_location: pathlib.Path = None) -> None:
    """Get a list of search tokens and find
    associated genomes. Perform fuzzy search
    on a search token that doesnt yield an
    exact match.

    Genomes with or without annotation  are
    written to the output file.

    Args:
        search_tokens (List[str]): A list of input search tokens.
        output_file (pathlib.Path): The output file. Can already exist and will be overwritten
        report_rich (bool, optional): If true, report also annotation. Only report genome names of false. Defaults to False.
    """    

    #### PARAMS ####
    evaluate_expression = False

    if db_location is None:
        db_location = mutils.DEFAULT_MOTUS_MGDB_LOCATION
    annodb_location = db_location / 'mOTUsv4.0.annotation.db'
    annodb_marker   = db_location / 'mOTUsv4.0.annotation.db.downloaded'
    mgdb_marker     = db_location / 'db_mOTU.downloaded'

    if not annodb_marker.exists():
        if not mgdb_marker.exists():
            logging.error('mOTUs marker gene database not downloaded. Download database with "motus downloadMGDB"')
            mutils.shutdown(1)
        logging.info('Need to download mOTUs annotation database (~17GB))')
        with urllib.request.urlopen(mutils.MOTUS_ANNODB_REMOTE_LOCATION) as response:
            total = int(response.info().get("Content-Length", -1))
            with open(str(annodb_location), "wb") as f, tqdm.tqdm(total=total, unit='B', unit_scale=True, desc='Downloading mOTUs annotation database') as pbar:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))

        DATABASE_PATH = str(annodb_location)
        annodb_marker.touch()
        logging.info('Finished downloading mOTUs annotation database.')
    else:
        DATABASE_PATH = str(annodb_location)
    #### END PARAMS ####

    search_tokens = sorted(set(search_tokens))
    search_token_2_genome_ids = {}
    resolver = SearchDB(DATABASE_PATH)
    for search_token in search_tokens:
        database_hits = resolver.exact_search_keyword(search_token)
        if len(database_hits) == 0:
            logging.info(f'Exact search for "{search_token}" didn\'t yield any hits.')
            logging.info('Trying with fuzzy search.')
            resolver.fuzzy_search_keyword(search_token)
            mutils.shutdown(1)
        else:
            logging.info(f'Exact search for "{search_token}" yielded {len(database_hits)} hit(s). Loading associated genomes.')
            genome_ids = resolver.get_genome_ids_by_database_hits(database_hits)
            logging.info(f'Found {len(genome_ids):,} associated genomes.')
            search_token_2_genome_ids[search_token] = genome_ids

    expression_2_genome_ids = {}
    if evaluate_expression:
        X = 0
        logging.info('Not implemented yet')
        exit(1)
        # here would be there expression eval
    else:
        expression_2_genome_ids = search_token_2_genome_ids
    if len(annotations_to_report) == 0:
        logging.info(f'Report mode: "basic". Writing names of genomes to {output_file}')
        with open(output_file, 'w') as handle:
            handle.write('GENOME\tQUERY\n')
            for expression, genome_ids in expression_2_genome_ids.items():
                genome_names = sorted(resolver.get_genome_names(genome_ids))
                for genome_name in genome_names:
                    handle.write(f'{genome_name}\t{expression}\n')
    else:
        logging.info(f'Collecting annotations {annotations_to_report} and writing to {output_file}')
        annotations_to_report = set(annotations_to_report)
        with open(output_file, 'w') as handle:
            handle.write(f'GENOME\tQUERY\t{'TAXONOMY\tmOTU\t' if 'TAXONOMY' in annotations_to_report else ''}{'EGGNOG\t' if 'EGGNOG' in annotations_to_report else ''}{'KEGG\t' if 'KEGG' in annotations_to_report else ''}\t{'PFAM\t' if 'PFAM' in annotations_to_report else ''}\n') #TODO add missing columns
            for expression, genome_ids in expression_2_genome_ids.items():
                logging.info(f'Search Token: {expression}:') 
                genomes = resolver.get_genomes_by_genome_ids(genome_ids)
                for genome in tqdm.tqdm(resolver.stringify_genomes(genomes), total=len(genomes), unit='genomes', desc='Writing genomes + annotations to file'):
                    
                    handle.write(f'{genome.get_genome_name()}\t{expression}\t{genome.get_tax_printable() + '\t' if 'TAXONOMY' in annotations_to_report else ''}{genome.get_motu_printable() + '\t' if 'TAXONOMY' in annotations_to_report else ''}{genome.get_eggnog_printable() + '\t' if 'EGGNOG' in annotations_to_report else ''}{genome.get_kegg_printable() + '\t' if 'KEGG' in annotations_to_report else ''}{genome.get_pfam_printable() + '\t' if 'PFAM' in annotations_to_report else ''}\n')

                
    resolver.close()
    


