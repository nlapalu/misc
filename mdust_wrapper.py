#!/usr/bin/env python

import subprocess
import tempfile
import sys
import os
import re
from optparse import OptionParser


class MdustWrapper(object):

    def __init__(self):
        self._options = None

    def stop_err(self, msg):
        sys.stderr.write("%s\n" % msg)
        sys.exit()

    def setAttributesFromCmdLine(self):
        description = "mdust_wrapper"
        description += "\nWrapper for mdust\n"
        description += "example: mdust_wrapper.py -i seq.fasta -v 27\n"
        parser = OptionParser(description = description, version = "0.1")
        parser.add_option("-i", "--input",	   dest = "FastaFile",   action = "store", type = "string", help = "Input Fasta File name [compulsory] [format: Fasta]", default = "")
        parser.add_option("-o", "--output",	dest = "outFile",   action = "store", type = "string", help = "output File name [compulsory] [format: fasta,tab or bed]", default = "")
        parser.add_option("-v", "--cutoff",	 dest = "cutoff",	 action = "store", type = "int",	help = "cutoff",	default = 28)
        parser.add_option("-w", "--wsize",	 dest = "wsize",	 action = "store", type = "int",	help = "window size",	default = 3)
        parser.add_option("-m", "--maskingletter",	 dest = "maskingletter",	 action = "store", type = "string",	help = "masking letter",	default = "N")
        parser.add_option("-f", "--format",	   dest = "format",  action = "store", type = "string", help = "format", default = "default")
        options = parser.parse_args()[0]
        self._setAttributesFromOptions(options)


    def _setAttributesFromOptions(self, options):
        self._options = options

        if self._options.FastaFile == "":
            raise Exception("Missing input file, please provide fasta file with -i file !")
        if self._options.outFile == "":
            raise Exception("Missing output file, please provide output file with -o file !")


    def run(self):

		prg = "mdust"
		args = ""
		args += " %s" % self._options.FastaFile
                args += " -v %d" % self._options.cutoff
                args += " -w %d" % self._options.wsize
                args += " -m %s" % self._options.maskingletter
                if self._options.format == "tab" or self._options.format == "bed":
                    args += " -c "
		cmd = "%s %s" %(prg, args)

		try:
			tmp_err = tempfile.NamedTemporaryFile().name
                        tmp_out = "outfile"
			tmp_stderr = open( tmp_err, 'wb' )
			tmp_stdout = open( tmp_out, 'wb' )
			proc = subprocess.Popen( args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr )
			returncode = proc.wait()
			tmp_stderr.close()
			# get stderr, allowing for case where it's very large
			tmp_stderr = open( tmp_err, 'rb' )
			tmp_stdout = open( tmp_out, 'rb' )

			stderr = ''
                        stdout = ''
			buffsize = 1048576
			try:
				while True:
					stdout += tmp_stdout.read( buffsize )
					if not stdout or len( stdout ) % buffsize != 0:
						break
			except OverflowError:
				pass
			tmp_stdout.close()

                        try:
				while True:
					stderr += tmp_stderr.read( buffsize )
					if not stderr or len( stderr ) % buffsize != 0:
						break
			except OverflowError:
				pass
			tmp_stderr.close()
			if stderr:
				raise Exception, stderr
		except Exception, e:
			self.stop_err( 'Error with mdust :\n' + str( e ) )

                if self._options.format == 'bed':
                    with open(tmp_out,"r") as fin:
                        with open(self._options.outFile, "w") as fout:
                            lineNumber = 0
                            for line in fin:
                                lineNumber += 1
                                m = re.search(r"^(\S+)\t(\d+)\t(\d+)\t(\d+)$", line)
                                if m is not None:
                                    fout.write("%s\t%d\t%d\n" % (m.group(1), int(m.group(3))-1, int(m.group(4))))
                                if m is None:
                                    raise Exception("\nLine %d '%s' does not has a mdust format." % (lineNumber, line))
                else:
                    os.rename(tmp_out,self._options.outFile)




if __name__ == "__main__":
	iWrapper = MdustWrapper()
	iWrapper.setAttributesFromCmdLine()
	iWrapper.run()
