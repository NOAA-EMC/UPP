.. attention:: 

   On Ursa, the usual process for generating a new/updated ``postxconfig*.txt`` file is slightly different due to a missing XML module (see `Issue #1250 <https://github.com/NOAA-EMC/UPP/issues/1250>`_). The following workaround has been developed:

   .. code-block:: console

      wget https://raw.githubusercontent.com/wiki/NOAA-EMC/UPP/perl_venv_create.sh
      chmod 755 perl_venv_create.sh
      ./perl_venv_create.sh perl_venv
      source perl_venv/bin/activate
      cpanm XML::LibXML
      cd /path/to/UPP
      cd parm
   
   Then, run ``PostXMLPreprocessor.pl`` as described above. 