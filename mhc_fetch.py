#!/usr/bin/env python

"""
Python script to fetch allele sequences from the IPD-MHC database

"""

import argparse
import requests
import os

API_ENDPOINT = 'https://www.ebi.ac.uk/cgi-bin/ipd/mhc/'

def fetchSequence(accession, type):

    data = requests.get(API_ENDPOINT + 'allele/' + accession + '/sequence?type=' + type)
    if data.status_code != 200:
        raise Exception('Unable to retrieve allele sequence!')

    data = data.json()

    return data['sequence']

if __name__ == "__main__":

    # Argument Handling
    ap = argparse.ArgumentParser(description=__doc__)

    # Input options
    io = ap.add_argument_group('Input options')
    io.add_argument('-g', '--gene', help='gene name', type=str, required=False)
    io.add_argument('-r', '--organism', help='organism name', type=str, required=False)

    # Output options
    oo = ap.add_argument_group('Output Options')
    oo.add_argument('-o', '--output', help='Output file.', type=str, required=False, default='./data')

    # Parse options
    cmd = ap.parse_args()

    # Generate the output folder
    if not os.path.exists(cmd.output):
        os.makedirs(cmd.output)

    # Fetch a list of alleles
    offset = 0;
    while True:

        # Prepare a list of arguments
        query = ''
        arguments = {'gene': cmd.gene, 'organism': cmd.organism, 'offset': str(offset), 'group': 'FISH'}
        for key, value in arguments.items():

            if value is None:
                continue

            query += '&' if query else '?'
            query += key + '=' + value

        # Send request
        data = requests.get(API_ENDPOINT + 'allele' + query)
        if data.status_code != 200:
            raise Exception('Unable to retrieve allele list!')

        # No more data to load
        data = data.json()
        if len(data['data']) == 0:
            break

        # Process each allele
        for value in data['data']:

            print('Fetching allele %s' % value['accession'])

            # Store allele infos
            with open(os.path.join(cmd.output, "IPD-MHC.txt"), 'w' if offset == 0 else 'a') as file:

                if offset == 0:
                    file.write('ACCESSION\tNAME\tSPECIES\n')

                file.write('{}\t{}\t{}\n'.format(value['accession'], value['name'], value['organism']))

            # Fetch nucleotide sequence
            with open(os.path.join(cmd.output, "IPD-MHC.nt"), 'w' if offset == 0 else 'a') as file:

                sequence = fetchSequence(value['accession'], 'full')

                # Append to file
                file.write('>MHC|{}|{}\n'.format(value['accession'], value['name'].replace('-','_')))
                file.write('{}\n'.format(sequence))
            

            # Fetch protein sequence
            with open(os.path.join(cmd.output, "IPD-MHC.aa"), 'w' if offset == 0 else 'a') as file:
            
                sequence = fetchSequence(value['accession'], 'protein')

                # Append to file
                file.write('>MHC|{}|{}\n'.format(value['accession'], value['name'].replace('-','_')))
                file.write('{}\n'.format(sequence))

            # Increase offset number
            offset += 1


    print('Successfully fetched %d alleles' % offset)
