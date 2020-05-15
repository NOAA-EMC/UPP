# This makefile is to create Post flat file from available and control XML file

FLATFILEDIR		= /glade/p/work/fossell/UPP/MERGE/DTC_TRUNK/parm

# Target flat file
FLATFILENAME		= postxconfig-NT.txt

# Source Post XML file
AVAILXMLFILENAME	= post_avblflds.xml
CTRLFILENAME		= postcntrl.xml

# Post flat file generator
PERLXML			= PostXMLPreprocessor.pl

# CMD with param
DTCPERLXML	 	= /usr/bin/perl $(PERLXML) $(CTRLFILENAME) $(AVAILXMLFILENAME) $(FLATFILENAME)

# File to look for change
XMLS = $(AVAILXMLFILENAME) $(CTRLFILENAME)

# If action is triggered; run the following
all: $(FLATFILENAME) 
$(FLATFILENAME): $(XMLS)
	$(DTCPERLXML)

clean:
	/bin/rm -f $(FLATFILENAME)
