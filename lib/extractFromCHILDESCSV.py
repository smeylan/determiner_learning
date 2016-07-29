import argparse
import os
import nltk
import itertools
import datetime
import codecs
import cPickle
import csv
import re

parser = argparse.ArgumentParser(description='Extract nouns and determiners from CSVs from CHILDES files')
parser.add_argument('-if','--inputfile', help='path of the file to read in', required=True)
parser.add_argument('-of','--outputfile', help='path of file to write', required=True)
parser.add_argument('-n','--noun', help='first or last noun flag', required=True)
parser.add_argument('-st','--stanfordTagger', help='path to the stanfor tagger directory', required=True) #'/Applications/stanford-postagger/'
args = parser.parse_args()

STANFORD_TAGGER_PKGDIR = args.stanfordTagger
STANFORD_TAGGER_MODEL = os.path.join(STANFORD_TAGGER_PKGDIR,'models/english-caseless-left3words-distsim.tagger')
STANFORD_TAGGER_JAR = os.path.join(STANFORD_TAGGER_PKGDIR,'stanford-postagger.jar')
NOUN_TAGS = set(['NN','NNP','NNS','NNPS'])
DET_WORDS = set(['a','an','the'])

if args.noun == 'FN':
    use_first_noun = True
elif args.noun == 'LN':
    use_first_noun = False
inputFile = args.inputfile
outputFile = args.outputfile


def extract_det_noun_phrases_v2(tagged,use_first_noun=True):
    """
    Extract ALL DET (non det, non-noun)* NOUN phrases from the tagged utt chunk, up to first or last NOUN

    Parameters:
    - tagged:
    - use_first_noun: if True (default), extracts from the determiner up to the FIRST noun after the determiner
                      if False, extracts up to the LAST noun after the determiner
    Returns:
    - phrases: The list of DET-NOUN phrases found
    - ranges: List of [start_tok_idx,end_tok_idx) tuples of chunk token indices, one tuple for each phrase found
    - danglers: List of booleans indicating whether the found phrases is "dangling" (ie. does not end with a NOUN)
    """
    toks, tags = zip(*tagged)
    phrases = list()
    ranges = list()
    danglers = list()
    i = 0
    while i < len(tags):
        if tags[i] == 'DT' and toks[i] in DET_WORDS: # we found the determiner, now iterate forward to the noun
            found_noun = False
            target_noun_idx = -1 # index of the noun to go with this determiner
            j = i+1
            while j < len(tags):
                if tags[j] in NOUN_TAGS: # we found a noun to go with the determiner
                    target_noun_idx = j
                    found_noun = True
                elif found_noun: # we have a noun and just encountered a non-noun, so end the (complete) phrase
                    break
                elif tags[j] == 'DT': # we don't have a noun, and we just encountered another determiner, so break
                    break
                j += 1
                if use_first_noun and found_noun: # If we're using the first noun, break here
                    break
            if found_noun:
                phrase = tagged[i:(target_noun_idx+1)]
                phrases.append(phrase)
                ranges.append((i,target_noun_idx+1))
                danglers.append(False)
            else:
                # print 'Dangling determiner {:s} for tagged chunk: {}'.format(det_word,tagged)
                dangle_phrase = tagged[i:(i+2)]
                phrases.append(dangle_phrase)
                ranges.append((i,(i+2)))
                danglers.append(True)
            #
            i = j
        else:
            i += 1
    assert len(phrases) == len(ranges)
    assert len(phrases) == len(danglers)
    return phrases, ranges, danglers


def export_pos_tagged_CHILDES(inputFile, outputFile, use_first_noun):
    """
    This is a test for Stephan - want to try the Stanford POS tagger and my phrase extraction code
    on the CHILDES data
    
    """ 
    
    # It's rather annoying that csvreader doesn't handle unicode, but they suggest this solution
    # when you have unicode data
    def unicode_csv_reader(unicode_csv_data, dialect=csv.excel, **kwargs):
        # csv.py doesn't do Unicode; encode temporarily as UTF-8:
        csv_reader = csv.reader(utf_8_encoder(unicode_csv_data),
                                dialect=dialect, **kwargs)
        for row in csv_reader:
            # decode UTF-8 back to Unicode, cell by cell:
            yield [unicode(cell, 'utf-8') for cell in row]
    
    
    def utf_8_encoder(unicode_csv_data):
        for line in unicode_csv_data:
            yield line.encode('utf-8')
    
    
    # The Penn Treebank tokenizer, which properly deals with contractions which is important
    # for the subsequent POS tagger
    toker = nltk.tokenize.TreebankWordTokenizer()
    
    # First load the glosses
    glosses = list()
    gloss_ids = list()
    with codecs.open(inputFile,'r','utf-8') as f:
        r = unicode_csv_reader(f)
        hdr = r.next()
        gloss_idx = hdr.index('Gloss')
        filename_idx = hdr.index('Filename')
        utt_num_idx = hdr.index('Utt.Number')
        spkr_idx = hdr.index('Speaker')
        
        assert gloss_idx >= 0 and utt_num_idx >= 0 and filename_idx >= 0
        prev_gloss_id = (None,None)
        unique_speakers = set()
        
        rows_loaded = 0
        print 'Filename: {:s}\n\tHeader: {:s}'.format(inputFile, hdr)
        for row in r:
            rows_loaded += 1
            # if rows_loaded > 1000:
            #     break
            
            chat_filename = row[filename_idx].strip().encode('utf-8') # just so we don't have u'xxx' for everything...
            utt_num = int(row[utt_num_idx].strip())
            gloss_id = (chat_filename,utt_num) # gloss id is really a combination of filename and utt number
            assert prev_gloss_id[0] != gloss_id[0] or prev_gloss_id[1] < gloss_id[1]                        
            
            prev_gloss_id = gloss_id
            gloss = row[gloss_idx].strip()
            speaker = row[spkr_idx].strip()
            unique_speakers.add(speaker)
                            
            if len(gloss) == 0:
                print 'Skipping zero length gloss: {}'.format(gloss_id)
                continue
            
            # Tokenize (using a tokenizer that's compatible with the POS tagger)
            gloss_toks = re.split('[ _]',gloss)
            gloss_toks = tuple(toker.tokenize(' '.join(gloss_toks))) # retokenize with the Penn Treebank tokenizer
            
            # Store the gloss info
            glosses.append(gloss_toks)
            gloss_ids.append(gloss_id)
            
        f.close()
        
        print 'Loaded {:d} rows, num glosses: {:d}'.format(rows_loaded,len(glosses))
        print 'Speaker column values: {}'.format(unique_speakers)
        
    
    # Next initialize and run the POS tagger on the glosses
    print 'Batch tagging glosses...'
    os.environ["CLASSPATH"] = STANFORD_TAGGER_PKGDIR+":"+STANFORD_TAGGER_PKGDIR+'/lib/*'
    os.environ["STANFORD_TAGGER_MODEL"] = STANFORD_TAGGER_MODEL
    pos_tagger = nltk.tag.StanfordPOSTagger(STANFORD_TAGGER_MODEL,STANFORD_TAGGER_JAR, encoding='utf-8')
    
    # Batch tagging...
    tagged_results = list()
    range_span = 1000
    for i in xrange(0,len(glosses),range_span):
        j = min(len(glosses),i+range_span)
        print 'Processing range [{:d},{:d}), {:.2f}% of total'.format(i,j,i*100.0 / len(glosses))
        partial_results = pos_tagger.tag_sents(glosses[i:j])
        tagged_results.extend(partial_results)
    
    if len(tagged_results) != len(glosses) or len(tagged_results) != len(gloss_ids):
        raise Exception('Lengths differ {:d} {:d} {:d}'.format(len(tagged_results),len(glosses),len(gloss_ids)))
        
    #     
    of = codecs.open(outputFile,'w','utf-8')
    for gloss_id,tagged,orig_gloss in itertools.izip(gloss_ids,tagged_results,glosses):
        unpacked_gloss,tags = zip(*tagged)        
        
        # Just make sure everything looks good here - that what we passed in to the POS tagger is what we get out
        if unpacked_gloss != orig_gloss:
            print 'orig_gloss: {} != unpacked tagged gloss: {}'.format(orig_gloss,unpacked_gloss)
        
        # Get the DET-NOUN phrases
        found_phrases = list()
        
        # Get the DET-NOUN phrases, version 1
        # if 'a' in orig_gloss:
        #     phrases,ranges,num_dangling = extract_det_noun_phrases(tagged,'a')
        #     found_phrases.extend(phrases)
        #     
        # if 'the' in orig_gloss:
        #     phrases,ranges,num_dangling = extract_det_noun_phrases(tagged,'the')
        #     found_phrases.extend(phrases)
        
        # Get the DET-NOUN phrases, version 2
        if 'a' in orig_gloss or 'the' in orig_gloss or 'an' in orig_gloss:
            phrases,ranges,dangle_bools = extract_det_noun_phrases_v2(tagged,use_first_noun=use_first_noun)
            phrases_good = [p for p,is_dangling in itertools.izip(phrases,dangle_bools) if not is_dangling]
            
            if len(phrases_good) == 0:
                # print "No phrases for input: {}, gloss id: {}".format(tagged,gloss_id)
                continue
            
            found_phrases.extend(phrases_good)
        
        for p in found_phrases:
            tagged_str = u' '.join([u'{:s}_{:s}'.format(x[0],x[1]) for x in p])
            gloss_str = u' '.join(orig_gloss)
            print >> of, u'{}\t{:s}\t{:s}'.format(gloss_id, gloss_str, tagged_str)
    
    of.close()
    print 'Done!'

def main():
    export_pos_tagged_CHILDES(inputFile, outputFile, use_first_noun)

if __name__ == '__main__':
    main()