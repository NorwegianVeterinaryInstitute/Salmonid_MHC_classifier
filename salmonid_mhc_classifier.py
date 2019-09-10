import os, sys, subprocess
from Bio import SeqIO
from ete3 import Tree, TreeStyle

def collect_database_info():
    table = {}
    for line in open('data/IPD-MHC.txt', 'r'):
        if 'ACCESSION' not in line:
            accession, name, species, previous_name = line.rstrip().split('\t')
            if len(name) == 11:
                proposed_name = name
            elif len(name) == 13:
                proposed_name = str(name)[0:11] + ':' + str(name)[11:13]
            elif len(name) == 15:
                proposed_name = name[0:11] + ':' + name[11:13] + ':' + name[13:15]
            table[accession] = [accession, name, species, previous_name, proposed_name]

    table_nt = {}
    for record in SeqIO.parse(open('data/IPD-MHC.nt', 'r'), 'fasta'):
        table_nt[record.id.split('|')[1]] = [record.id, record.name, record.description, str(record.seq)]

    table_aa = {}
    for record in SeqIO.parse(open('data/IPD-MHC.aa', 'r'), 'fasta'):
        table_aa[record.id.split('|')[1]] = [record.id, record.name, record.description, str(record.seq)]
        
    print('Database information collected')
    print()
    return(table, table_nt, table_aa)
    
def execute(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    p_status = p.wait()

def check_nt_length(seq_len):
    if seq_len % 3 != 0:
        print (name + ' length is not divisble by 3. Check!!')
    else:
        return('ok')
    
def write_nt_file(record, nt_file):
    out_file = open(nt_file, 'w')
    SeqIO.write(record, out_file, 'fasta')
    out_file.close()

def translate(name, nt_file, seq_len, aa_file):
    execute('transeq -sequence ' + nt_file + ' -outseq ' + aa_file)
    for fasta in SeqIO.parse(open(aa_file, 'r'), 'fasta'):
        aa_len = len(str(fasta.seq))
        if (seq_len / 3) !=  aa_len:
            print (name + 'is not completed translated to protein sequence. Check!!')
        else:
            return('ok')
        
def muscle_aln(name, nt_file, input_file, specific_db, aln_file, tree_file, tree_png_file):
    execute('cat ' + nt_file + ' data/' + specific_db + '.nt > ' + input_file)
    execute('muscle -in ' + input_file + ' -out ' + aln_file + ' -tree2 ' + tree_file)
    Tree(tree_file).render(tree_png_file, dpi=200, w=400, units='mm')
    print(name + ' (2/5 tasks): Muscle complete')
    
def clade_extraction(name, tree_file):
    siblings = []
    node = Tree(tree_file).search_nodes(name=name)[0]
    for i in range(4):
        for line in str(node).split():
            potential_sibling_id = line.split('-')[-1]
            if len(potential_sibling_id) >= 2 and name not in potential_sibling_id:
                siblings.append(potential_sibling_id.split('|')[1])
        node = node.up
    return(set(siblings))


def clade_seq_extraction(clade_id, clade_xx_file, table_xx):
    out_file = open(clade_xx_file , 'w')
    out_file.write('>' + clade_id + '\n')
    out_file.write(table_xx[clade_id][3] + '\n')
    out_file.close()

def water_aln(name, xx_file, clade_xx_file, water_xx_out_file):
    execute('water -asequence ' + xx_file + ' -bsequence ' + clade_xx_file + ' -gapopen 10 -gapextend 0.5 -outfile ' + water_xx_out_file)

    print(name + ' (4-5/5 tasks): Water comparison complete')
    
    for line in open(water_xx_out_file, 'r'):
        if 'Identity' in line:
            identity = line.replace('# ', '').rstrip()
        if 'Similarity' in line:
            similarity = line.replace('# ', '').rstrip()

    return(identity, similarity)
 
def main(input_file, report_file_name, specific_db):
    (table, table_nt, table_aa) = collect_database_info()
    
    report_file = open(report_file_name, 'w')
    
    execute('rm -r scratch output')
    
    execute('mkdir scratch')
    execute('mkdir output')
    
    for record in SeqIO.parse(open(input_file, 'r'), 'fasta'):
        name = record.id
        seq = str(record.seq)
        seq_len = len(seq) 
        
        nt_file = 'output/' + name + '.fa'
        aa_file = 'output/' + name + '.aa'
        
        input_file = 'scratch/' + name + '_muscle_input.nt'
        aln_file = 'output/' + name + '_muscle_output.aln.txt'
        tree_file = 'output/' + name + '_muscle_output.tree.txt'
        tree_png_file = 'output/' + name + '_muscle_output.tree.png'
#        tree_svg_file = 'scratch/' + name + '_muscle_output.tree.svg'
#        tree_pdf_file = 'output/' + name + '_muscle_output.tree.pdf'
                
        if check_nt_length(seq_len) is not 'ok':
            break
           
        write_nt_file(record, nt_file)
        
        if translate(name, nt_file, seq_len, aa_file) is not 'ok':
            break
        else:
            print(name + ' (1/5 tasks): Translate complete')
        
        muscle_aln(name, nt_file, input_file, specific_db, aln_file, tree_file, tree_png_file)
        
#        clade_id = clade_extraction(name, tree_file)
        clade_ids = clade_extraction(name, tree_file)
        print (clade_ids)
        for clade_id in clade_ids:
            clade_nt_file = 'scratch/' + clade_id + '.nt'
            clade_aa_file = 'scratch/' + clade_id + '.aa'
        
            water_nt_out_file =  'output/' + name + '_' + clade_id + '_water_nt.txt'
            water_aa_out_file =  'output/' + name + '_' + clade_id + '_water_aa.txt'
        
            clade_seq_extraction(clade_id, clade_nt_file, table_nt)
            clade_seq_extraction(clade_id, clade_aa_file, table_aa)
        
            (nt_identity, nt_similarity) = water_aln(name, nt_file, clade_nt_file, water_nt_out_file)
            (aa_identity, aa_similarity) = water_aln(name, aa_file, clade_aa_file, water_aa_out_file)
            
            report_file.write('*********************************************************************************\n')
        
            report_file.write(name + '\t' + str(seq_len) + '\t'+ str(seq_len/3) + '\n')
            report_file.write(clade_id + '\n')  
            report_file.write('\t'.join(table[clade_id]) + '\n')
            report_file.write('nt: ' + nt_identity + '\t' + nt_similarity + '\n')
            report_file.write('aa: ' + aa_identity + '\t' + aa_similarity + '\n')
        
        report_file.write('*********************************************************************************\n')
        report_file.write('\n')
        print()
    
    report_file.write('*********************************************************************************\n')
    report_file.close()

if __name__ == '__main__':
    input_file = sys.argv[1]
    specific_db = sys.argv[2]
    report_file_name = sys.argv[3]
    main(input_file, report_file_name, specific_db)

    print('All done')
