#!/usr/bin/env python

import sys
import os
import logging
import argparse
import pysam

class ExtractSeqFromVCF(object):

    def __init__(self, logLevel='ERROR'):
        """__init__"""

        self.dSamples = {}
        self.logLevel = logLevel
        logging.basicConfig(level=self.logLevel)
        
    def __del__(self):
        pass

    def readListOfVarFiles(self, lFiles):
        """Read list of Variant Files"""

        lVarFiles = []
        try:
            with open(lFiles, 'r') as f:
                for line in f:
                    if line:
                        lVarFiles.append(line.rstrip())
            return lVarFiles

        except(Exception):
            logging.error('Error when reading input file list')
            sys.exit(1)

    def isVariantFileIndexed(self, varFile):
        """Check existing index"""

        if os.path.isfile('{}.tbi'.format(varFile)) or os.path.isfile('{}.csi'.format(varFile)):
            return True
        else:
            logging.info('No index found for variant file: {}'.format(varFile))
            return False
        
    def indexingVariantFile(self, varFile):
        """Index variant file with Tabix"""

        logging.info('Trying to index file: {}'.format(varFile))
        pysam.tabix_index(varFile,force=True,preset="vcf")
        if not self.isVariantFileIndexed(varFile):
            raise Exception("Can not index file: {}".format(varFile))
        return True

    def addSamples(self, lSamples, varFile):
        """Add samples to list of Samples, check unicity"""

        logging.debug('Adding new samples')
        for id in lSamples:
            if id not in self.dSamples:
                self.dSamples[id] = Sample(id)
            else:
                logging.error('Sample: {} is already in sample list, please change the name in the variant file: {}'.format(id, varFile))
                sys.exit(1) 
        return lSamples 
        
    def checkIndexOfVariantFiles(self, lVarFiles):
        """Check index for all variant files"""

        try:
            for varFile in lVarFiles:
               if not self.isVariantFileIndexed(varFile):
                   self.indexingVariantFile(varFile)
        except Exception as e:
            logging.error(e.message)
            logging.info('Missing variant file index, extraction will not be possible')
            sys.exit(1)

    def extract(self, contig, start, end, lVarFiles):
        """Extract region with 0-based indexing system"""
        
        logging.info('Extracting from variant files')

        self.checkIndexOfVariantFiles(lVarFiles)
         
        for varFile in lVarFiles:
            varFileHeader = pysam.VariantHeader()
            f = pysam.VariantFile(varFile, header=varFileHeader)

            # get list of samples
            lCurrentSamples = self.addSamples(list(f.header.samples), varFile)

            # get list of pos
            for val in f.fetch(str(contig), int(start), int(end)):
               for sample in val.samples.items():
                   lAlleles, depth = self.getAllelesFromPosition(list(val.alleles), dict(val.info), sample[1].items())
                   if depth > 0:
                       genotype = (self._getFormatTagValue(sample[1].items(),'GT',index=0),self._getFormatTagValue(sample[1].items(),'GT',index=1))
                   else:
                       genotype = (None, None)
                  # depth = self._getFormatTagValue(sample[1].items(), 'DP')
                   self.dSamples[sample[0]].addPosition(Position(val.contig, val.pos, genotype=genotype, depth=depth, lAlleles=lAlleles, lFilters=val.filter))

    def getAllelesFromPosition(self, lbases, dPosInfoTags, lPosFormatTags):
        """get Allele from position"""

        lAlleles = []
        baseRefAllele = lbases.pop(0)

        # fix bug variant with no mapping
        depthTotal = self._getInfoTagValue(dPosInfoTags, 'DP')
        if depthTotal == 0:
            lAlleles.append(Allele(baseRefAllele, 0, 0, refAllele=True))
            lAlleles.append(Allele(".", 0, 0, refAllele=False))
            return lAlleles, depthTotal

        #ref
        depthAllAlleles = self._getFormatTagValue(lPosFormatTags, 'DP')
        depthRefAllele = self._getFormatTagValue(lPosFormatTags, 'RO')
        freqRefAllele =  float(depthRefAllele) / depthAllAlleles
        lAlleles.append(Allele(baseRefAllele,freqRefAllele,depthRefAllele,refAllele=True))
        #alts
        for i,base in enumerate(lbases):
            depthAltAllele = self._getFormatTagValue(lPosFormatTags, 'AO', index=i)
            baseAltAllele = base
            freqAltAllele =  float(depthAltAllele) / depthAllAlleles
            lAlleles.append(Allele(baseAltAllele,freqAltAllele,depthAltAllele,refAllele=False))

        return lAlleles, depthTotal

    def _getFormatTagValue(self, lTags, tag, index=0):
        """return specific format tag value, index specific"""

        for t in lTags:
            if t[0] == tag:
                if type(t[1]) == tuple:
                    return t[1][index]
                else:
                    return t[1]

    def _getInfoTagValue(self, dTags, tag, index=0):
        """return specific info tag value, index specific"""

        if type(dTags[tag]) == tuple:
            return dTags[tag][index]
        else:
            return dTags[tag]

    def export(self, format, contig, start, end, freq, nbAlleles, minDP, maxDP, exportFile):
        """export results"""

        if format == 'fasta':
            self.exportToFasta(contig, start, end, freq, nbAlleles, minDP, maxDP, exportFile)
        elif format == 'algmt':
            self.exportToAlgmt(contig, start, end, freq, nbAlleles, minDP, maxDP, exportFile)


    def exportToFasta(self, contig, start, end, freq, nbAlleles, minDP, maxDP, exportFile=None, charLen=60):
        """Export to Fasta format"""

        dExportSamples = {}

        for sample in self.dSamples.keys():
            index = 0
            allPos = range(start, end+1)
            dPos = dict(zip(allPos, ['N']*len(allPos)))
            dRefPos = dict(zip(allPos, ['N']*len(allPos)))

            for i,pos in enumerate(range(start, end+1)):
                if self.dSamples[sample].getPosition(contig, pos):
                    if not self.isPositionFiltered(self.dSamples[sample].getPosition(contig, pos), freq, nbAlleles, minDP, maxDP):
                        allele = self.dSamples[sample].getPosition(contig,pos).getMajorAllele().base
                        ref = self.dSamples[sample].getPosition(contig,pos).getRefAllele().base
                        maxNbBases = max(len(allele),len(ref))
                        lbases = ['-'] * maxNbBases
                        lrefbases = ['-'] * maxNbBases
                        lbases[0:len(allele)] = list(allele)
                        lrefbases[0:len(ref)] = list(ref)
                        for j,base in enumerate(lrefbases):
                            if base != '-':
                                dRefPos[pos+j] = lrefbases[j]
                                dPos[pos+j] = lbases[j]
                            else:
                                dPos[pos+j-1] = ''.join(lbases[j-1:])
                                break
            for pos in dPos:
                dPos[pos] = [x for x in dPos[pos] if x != '-']
            dExportSamples[sample] = dPos

        if exportFile:
            try:
                with open(exportFile, 'w') as f:
                    for sample in dExportSamples:
                        f.write(self._formatToFasta('{}:{}:{}:{}'.format(sample, contig, start,end), [dExportSamples[sample][pos] for pos in range(start, end+1)], 60))
            except Exception as e:
                logging.error(e)
        else:
            print("\n### sequence in fasta format ###\n") 
            for sample in dExportSamples:
                print(self._formatToFasta('{}:{}:{}:{}'.format(sample, contig, start,end), [dExportSamples[sample][pos] for pos in range(start, end+1)], 60))
 

    def exportToAlgmt(self, contig, start, end, freq, nbAlleles, minDP, maxDP, exportFile=None):
        """Export to Algmt format"""

        dExportSamples = {}
        dExportRefs = {}

        for sample in self.dSamples.keys():

            allPos = range(start, end+1)
            dPos = dict(zip(allPos, ['N']*len(allPos)))
            dRefPos = dict(zip(allPos, ['N']*len(allPos)))

            for i,pos in enumerate(range(start, end+1)):
                if self.dSamples[sample].getPosition(contig, pos):
                    if not self.isPositionFiltered(self.dSamples[sample].getPosition(contig, pos), freq, nbAlleles, minDP, maxDP):
                        allele = self.dSamples[sample].getPosition(contig,pos).getMajorAllele().base
                        ref = self.dSamples[sample].getPosition(contig,pos).getRefAllele().base
                        maxNbBases = max(len(allele),len(ref))
                        lbases = ['-'] * maxNbBases
                        lrefbases = ['-'] * maxNbBases
                        lbases[0:len(allele)] = list(allele)
                        lrefbases[0:len(ref)] = list(ref)
                        for j,base in enumerate(lrefbases):
                            if base != '-':
                                dRefPos[pos+j] = lrefbases[j]
                                dPos[pos+j] = lbases[j]
                            else:
                                dPos[pos+j-1] = ''.join(lbases[j-1:])
                                break

            dExportSamples[sample] = dPos
            dExportRefs[sample] = dRefPos


        dConsRefSeq = self._getConsensusRefSeq(dExportRefs, start, end)
        
        dAlignedRef = {}
        dAlignedSamples = dict(zip(self.dSamples.keys(),['']*len(self.dSamples.keys())))
        for sample in dAlignedSamples:
            dAlignedSamples[sample] = {}

        for pos in range(start, end+1):
            maxNbBases = max([len(dExportSamples[sample][pos]) for sample in self.dSamples.keys()]) 
            lrefbases = ['-'] * maxNbBases
            lrefbases[0:len(dConsRefSeq[pos])] = list(dConsRefSeq[pos])
            dAlignedRef[pos] = ''.join(lrefbases)
            for sample in dExportSamples:
                lbases = ['-'] * maxNbBases
                lbases[0:len(dExportSamples[sample][pos])] = list(dExportSamples[sample][pos])
                dAlignedSamples[sample][pos] = ''.join(lbases)

        if exportFile:
            try:
                with open(exportFile, 'w') as f:
                    f.write(self._formatToFasta('reference:{}:{}:{}'.format(contig, start,end), [dAlignedRef[pos] for pos in range(start, end+1)], 60)) 

                    for sample in dAlignedSamples:
                        f.write(self._formatToFasta('{}:{}:{}:{}'.format(sample, contig, start,end), [dAlignedSamples[sample][pos] for pos in range(start, end+1)], 60))
    
            except Exception as e:
                logging.error(e)
        else:
            print("\n### alignment in fasta format ###\n") 
            print(self._formatToFasta('reference:{}:{}:{}'.format(contig, start,end), [dAlignedRef[pos] for pos in range(start, end+1)], 60))

            for sample in dAlignedSamples:
                print(self._formatToFasta('{}:{}:{}:{}'.format(sample, contig, start,end), [dAlignedSamples[sample][pos] for pos in range(start, end+1)], 60))
    
    def _formatToFasta(self, header, seq, charLen):
        """ print fasta """

        fseq = '>{}\n'.format(header)
        l = []
        idx = 0
        for i in seq:
            l.extend(list(i))
            idx += len(i)
            if idx >= charLen:
                fseq += '{}\n'.format(''.join(l[0:charLen]))
                l = l[charLen:]
                idx = idx - charLen
        if l:
            fseq += '{}\n'.format(''.join(l))

        return fseq 
            
    def _getConsensusRefSeq(self, dExportRefs, start, end):
        """..."""

        allPos = range(start, end+1)
        dRefPos = dict(zip(allPos, ['N']*len(allPos)))

        for pos in range(start, end+1):
            base = 'N'
            for sample in dExportRefs:
                if dExportRefs[sample][pos] != 'N':
                    if base == 'N':
                        base = dExportRefs[sample][pos]
                    elif base != dExportRefs[sample][pos] and dExportRefs[sample][pos] != 'N':
                        raise Exception('Could not get a reliable consensus \
                        for reference sequence; position {} get {} or {} \
                        base'.format(base, dExportRefs[sample][pos]))
                dRefPos[pos] = base
                
        return dRefPos    
        
                
    def isPositionFiltered(self, pos, freq, nbAlleles, minDP, maxDP):
        """Test if position must be filtered"""

        if len(pos.lAlleles) > nbAlleles:
            return True
        if pos.getMajorAllele().freq < freq:
            return True
        if minDP and pos.depth < minDP:
            return True
        if maxDP and pos.depth > maxDP:
            return True
        return False 
                
    def exportToAlgmtOld(self, contig, start, end):
        """Export to Algmt format"""

        lSampleNames = []
        for sample in sorted(self.dSamples.keys()):
            lSampleNames.append(sample)

        lSampleBases = []
        for i in range(0, len(lSampleNames)):
            lSampleBases.append([])

        lRefBases = ['.'] * (end+1 - start)  

        for pos in range(start, end+1):
            for i,sample in enumerate(sorted(self.dSamples.keys())):
                if self.dSamples[sample].getPosition(contig, pos):
                    lSampleBases[i].append(self.dSamples[sample].getPosition(contig,pos).getMajorAllele().base)
                    refBases = self.dSamples[sample].getPosition(contig,pos).getRefAllele().base
                    lRefBases[i:len(refBases)] = list(refBases)
                else:
                    lSampleBases[i].append(".")
                    lRefBases[i].append(".")

        #print lSampleBases
        lRefBases = []

        for pos in range(start, end+1):
            lRefBases.append('*')

        lDisplayRefBases = ['']
        lDisplaySampleBases = ['']*len(lSampleNames)
        #print lDisplayBases
        size = 60
        index = 0
        for i,base in enumerate(lRefBases):
            #for sample in lSampleBases:
            #  print sample
            nbBaseMax = max([len(sample[i]) for sample in lSampleBases])
            #print nbBaseMax
            for j,sample in enumerate(lSampleBases):
                bases = ['-'] * nbBaseMax 
                bases[0:len(sample[i])] = list(sample[i])
                lDisplayBases[j] += "".join(bases)          
                
        for sample in lDisplayBases:
            print sample        
            


class Allele(object):

    def __init__(self, base, freq, depth, refAllele=False):

        self.base = base
        self.freq = freq
        self.depth = depth
        self.isRefAllele = refAllele
        

    def __gt__(self, other):
        """Comparison greater than"""

        return self.freq > other.freq 

    def __repr__(self):
        """__repr__"""

        return "base:{},freq:{},depth:{},ref:{}".format(self.base,self.freq, self.depth, self.isRefAllele)


class Position(object):

    def __init__(self, ref, pos, genotype=None, depth=None, lFilters=[], lAlleles=[]):
        """__init__"""
        
        self.ref = ref
        self.pos = pos
        self.genotype = genotype
        self.depth = depth
        self.lFilters = lFilters
        self.lAlleles = lAlleles

    def getMajorAllele(self):
        """return Major Allele"""

        return max(self.lAlleles)

    def getRefAllele(self):
        """return Reference Allele"""

        lAlleles = [allele for allele in self.lAlleles if allele.isRefAllele == True]
        if lAlleles:
            return lAlleles[0]
        else:
            return None


class Sample(object):

    def __init__(self, id):
        """__init__"""

        self.id = id
        self.lPositions = []
        self.dIndexPositions = {}


    def addPosition(self, pos):
        """add a position to this sample"""

        self.lPositions.append(pos)
        self.dIndexPositions[(pos.ref,pos.pos)] = len(self.lPositions) - 1
        return len(self.lPositions)

    def getPosition(self, ref, pos):
        """return a position if in lPositions"""

        if (ref, pos) in self.dIndexPositions:
            return self.lPositions[self.dIndexPositions[(ref, pos)]]
        else:
            return None


if __name__ == '__main__':

    program = 'get_sequence_from_vcf'
    version = 0.1
    description = 'Extract most probable sequence from vcf entries. Only major \
                   allele is considered, No phasing ! This tool is reliable \
                   for haploid genomes, not for polyploid ones.'
    parser = argparse.ArgumentParser(prog=program)
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("contig", help="contig", type=str)
    parser.add_argument("start", help="Start coordinate, 1-based indexing system", type=int)
    parser.add_argument("end", help="End coordinate, 1-based indexing system", type=int)
    parser.add_argument("listFH", help="file with list of vcf files to analyze", type=str)

    parser.add_argument("-f","--freq", help="minimal frequence to consider the allele, filtered under this value [default=0.0]", type=float, default=0.0)
    parser.add_argument("-n","--nbAlleles", help="maximal number of alleles allowed per sample, filtered upper this value [default=2]", type=int, default=2)
    parser.add_argument("-mi","--minDP", help="minimal depth of coverage to consider the allele, filtered under this value, [default=1]", type=int, default=1)
    parser.add_argument("-ma","--maxDP", help="maximal depth of coverage to consider the allele, filtered upper this value, [default=0]", type=int, default=0)
    parser.add_argument("-sf","--seqFasta", help="export sequence in fasta file", type=str, default=None)
    parser.add_argument("-af","--algmtFasta", help="export alignment in fasta file", type=str, default=None)
    parser.add_argument("--version", action='version', version='{} {}'.format(program,version))
    parser.add_argument("-v", "--verbosity", type=int, choices=[1,2,3],
                        help="increase output verbosity 1=error, 2=info, 3=debug")

    args = parser.parse_args()

    logLevel='ERROR'
    if args.verbosity == 1:
        logLevel = 'ERROR'
    if args.verbosity == 2:
        logLevel = 'INFO'
    if args.verbosity == 3:
        logLevel = 'DEBUG'
    logging.getLogger().setLevel(logLevel)

    try:
        if args.freq < 0 or args.freq > 1:
            raise ValueError("Error for frequency value, expected value: 0 < x < 1, ex 0.7") 

        if args.nbAlleles < 1:
            raise ValueError("Error, nbAlleles must be >= 1")

        if args.minDP < 0 or args.maxDP < 0 :
            raise ValueError("Error, minDP and maxDP must be > 0")
        if args.minDP != 0 and args.maxDP != 0 :
            if args.minDP > args.maxDP:
                raise ValueError ("Error, minDP and maxDP must be > 0")

    except ValueError as e:
        logging.error('ValueError Exception: {}'.format(e.message))
        sys.exit(1)
        

    i = ExtractSeqFromVCF(logLevel)
    lVarFiles = i.readListOfVarFiles(args.listFH)
    i.extract(args.contig, args.start-1, args.end, lVarFiles)
    i.export('fasta', args.contig, args.start, args.end, args.freq, args.nbAlleles, args.minDP, args.maxDP, args.seqFasta)    
    i.export('algmt', args.contig, args.start, args.end, args.freq, args.nbAlleles, args.minDP, args.maxDP, args.algmtFasta) 
