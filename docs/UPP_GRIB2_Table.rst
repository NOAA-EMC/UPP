================================
GRIB2 Fields Produced by Unipost
================================

GRIB2 fields produced by *unipost* (column 1), abbreviated names
used in the *postcntrl.xml* file (column 2), corresponding standard
grib2 pname (column 3), corresponding grib identification number for the
vertical coordinate (column 4), and corresponding array location UPP
uses to store the variable in parallel arrays (column 5).

.. csv-table::
   :file: /d1/hertneky/upp/docs/UPP_GRIB2_Table.csv
   :header: "Field Description","Name in Grib2 Control File","Grib2 pname","Vertical Level","UPP ID"
   :widths: 35,35,10,10,10

\*See Appendix A of the UPP Users Guid

\***4 types of CAPE and CIN can be output with use of the Levels control
line in the wrf\_cntrl.parm (nmb\_cntrl.parm file).

Surface based CAPE/CIN is output at one grib record, while the remaining
three types are output within one grib record in 3 levels.

    Level 1: Surface Based CAPE/CIN

    Level 2: Best Boundary Layer CAPE/CIN

    Level 3: Mixed Layer CAPE/CIN

    Level 4: Most Unstable CAPE/CIN
