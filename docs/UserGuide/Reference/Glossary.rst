.. _Glossary:

**********
Glossary
**********

.. glossary::

   CAPE
      Convective Available Potential Energy. 

   CCPP
      The `Common Community Physics Package <https://dtcenter.org/community-code/common-community-physics-package-ccpp>`_ is a forecast-model agnostic, vetted collection of code containing atmospheric physical parameterizations and suites of parameterizations for use in Numerical Weather Prediction (:term:`NWP`) along with a framework that connects the physics to the host forecast model.

   CIN
      Convective Inhibition.

   CRTM
      The `Community Radiative Transfer Model <https://www.jcsda.org/jcsda-project-community-radiative-transfer-model>`__ (CRTM) is a fast and accurate radiative transfer model developed at the `Joint Center for Satellite Data Assimilation <https://www.jcsda.org/>`__ (JCSDA) in the United States. It is a sensor-based radiative transfer model and supports more than 100 sensors, including sensors on most meteorological satellites and some from other remote sensing satellites. 

   Component
      A software element that has a clear function and interface. In Earth system models, components are often single portions of the Earth system (e.g., atmosphere, ocean, or land surface) that are assembled to form a whole.

   Component Repository
      A :term:`repository` that contains, at a minimum, source code for a single component.

   CONUS
      Continental United States

   CAM
   convection-allowing models
      Convection-allowing models (CAMs) are models that run on high-resolution grids (usually with grid spacing at 4km or less). They are able to resolve the effects of small-scale convective processes. They typically run several times a day to provide frequent forecasts (e.g., hourly or subhourly). 

   data assimilation
      Data assimilation is the process of combining observations, model data, and error statistics to achieve the best estimate of the state of a system. One of the major sources of error in weather and climate forecasts is uncertainty related to the initial conditions that are used to generate future predictions. Even the most precise instruments have a small range of unavoidable measurement error, which means that tiny measurement errors (e.g., related to atmospheric conditions and instrument location) can compound over time. These small differences result in very similar forecasts in the short term (i.e., minutes, hours), but they cause widely divergent forecasts in the long term. Errors in weather and climate forecasts can also arise because models are imperfect representations of reality. Data assimilation systems seek to mitigate these problems by combining the most timely observational data with a "first guess" of the atmospheric state (usually a previous forecast) and other sources of data to provide a "best guess" analysis of the atmospheric state to start a weather or climate simulation. When combined with an "ensemble" of model runs (many forecasts with slightly different conditions), data assimilation helps predict a range of possible atmospheric states, giving an overall measure of uncertainty in a given forecast.

   dycore
   dynamical core
      Global atmospheric model based on fluid dynamics principles, including Euler's equations of motion.

   echo top
      The radar-indicated top of an area of precipitation. Specifically, it contains the height of the 18 dBZ reflectivity value.

   EMC
      The `Environmental Modeling Center <https://www.emc.ncep.noaa.gov/emc_new.php>`__. 

   EPIC
      The `Earth Prediction Innovation Center <https://epic.noaa.gov/>`__ seeks to accelerate scientific research and modeling contributions through continuous and sustained community engagement in order to produce the most accurate and reliable operational modeling system in the world. 

   ESG
      Extended Schmidt Gnomonic (ESG) grid. The ESG grid uses the map projection developed by Jim Purser of NOAA :term:`EMC`. 

   ESMF
      `Earth System Modeling Framework <https://earthsystemmodeling.org/docs/release/latest/ESMF_usrdoc/>`__. The ESMF defines itself as “a suite of software tools for developing high-performance, multi-component Earth science modeling applications.” 

   FV3
      The Finite-Volume Cubed-Sphere :term:`dynamical core` (dycore). Developed at NOAA's `Geophysical 
      Fluid Dynamics Laboratory <https://www.gfdl.noaa.gov/>`__ (GFDL), it is a scalable and flexible dycore capable of both hydrostatic and non-hydrostatic atmospheric simulations. It is the dycore used in the UFS Weather Model.

   GFS
      The `Global Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast>`__. The GFS is a National Centers for Environmental Prediction (:term:`NCEP`) weather forecast model that generates data for dozens of atmospheric and land-soil variables, including temperatures, winds, precipitation, soil moisture, and atmospheric ozone concentration. The system couples four separate models (atmosphere, ocean, land/soil, and sea ice) that work together to accurately depict weather conditions.

   GRIB2 
      The second version of the World Meterological Organization's (WMO) standard for distributing gridded data.  

   GSI
      `Gridpoint Statistical Interpolation <https://dtcenter.org/community-code/gridpoint-statistical-interpolation-gsi>`__ (GSI) is a variational data assimilation system, designed to be flexible, state-of-art, and run efficiently on various parallel computing platforms. It supports :term:`RRFS` features. GSI code is publicly available `on GitHub <https://github.com/NOAA-EMC/GSI>`__, and fix file data is publicly available `here <https://ftp.emc.ncep.noaa.gov/jcsda/WDQMS/NCEP/GSI-FIX/>`__. 

   HPC-Stack
      `HPC-Stack <https://github.com/JCSDA/spack-stack>`__ is a repository that provides a unified, shell script-based build system for building the software stack required for numerical weather prediction (NWP) tools such as the `Unified Forecast System (UFS) <https://ufscommunity.org/>`__ and the `Joint Effort for Data assimilation Integration (JEDI) <https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/>`__ framework. It is being phased out in favor of :term:`spack-stack`. `HPC-Stack documentation <https://hpc-stack-epic.readthedocs.io/en/latest/>`__ is available, but the repository and documentation is rarely updated since it is being deprecated.

   HRRR
      `High Resolution Rapid Refresh <https://rapidrefresh.noaa.gov/hrrr/>`__. The HRRR is a NOAA real-time 3-km resolution, hourly updated, cloud-resolving, convection-allowing atmospheric model initialized by 3-km grids with 3-km radar assimilation. Radar data is assimilated in the HRRR every 15 min over a 1-hour period adding further detail to that provided by the hourly data assimilation from the 13-km radar-enhanced Rapid Refresh.

   JCSDA
   Joint Center for Data Satellite Assimilation
      The `Joint Center for Satellite Data Assimilation <https://www.jcsda.org/>`__ is a multi-agency research center hosted by the University Corporation for Atmospheric Research (`UCAR <https://www.ucar.edu/>`__). JCSDA is dedicated to improving and accelerating the quantitative use of research and operational satellite data in weather, ocean, climate, and environmental analysis and prediction systems.

   LAM
      Limited Area Model (grid type), formerly known as the "Stand-Alone Regional" or SAR. LAM grids use a regional (rather than global) configuration of the :term:`FV3` :term:`dynamical core`. 

   MPI
      MPI stands for Message Passing Interface. An MPI is a standardized communication system used in parallel programming. It establishes portable and efficient syntax for the exchange of messages and data between multiple processors that are used by a single computer program. An MPI is required for high-performance computing (HPC) systems.

   MRW
   Medium-Range Weather Application
      The `Medium-Range Weather Application <https://github.com/ufs-community/ufs-mrweather-app>`__ is a UFS Application that targets predictions of atmospheric behavior out to about two weeks. It packages a prognostic atmospheric model (the UFS Weather Model), pre- and post-processing tools, and a community workflow.

   NAM
      `North American Mesoscale Forecast System <https://www.ncei.noaa.gov/products/weather-climate-models/north-american-mesoscale>`_. NAM generates multiple grids (or domains) of weather forecasts over the North American continent at various horizontal resolutions. Each grid contains data for dozens of weather parameters, including temperature, precipitation, lightning, and turbulent kinetic energy. NAM uses additional numerical weather models to generate high-resolution forecasts over fixed regions, and occasionally to follow significant weather events like hurricanes.

   namelist
      A namelist defines a group of variables or arrays. Namelists are an I/O feature for format-free input and output of variables by key-value assignments in Fortran compilers. Fortran variables can be read from and written to plain-text files in a standardised format, usually with a ``.nml`` file ending.

   NCAR
      The `National Center for Atmospheric Research <https://ncar.ucar.edu/>`__. 

   NCEP
      National Centers for Environmental Prediction (NCEP) is an arm of the National Weather Service
      consisting of nine centers. More information can be found at https://www.ncep.noaa.gov.

   NEMSIO
      A binary format for atmospheric model output from :term:`NCEP`'s Global Forecast System (:term:`GFS`).

   netCDF
      NetCDF (`Network Common Data Form <https://www.unidata.ucar.edu/software/netcdf/>`__) is a file format and community standard for storing multidimensional scientific data. It includes a set of software libraries and machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data.

   NUOPC
      The `National Unified Operational Prediction Capability <https://earthsystemmodeling.org/nuopc/>`__ Layer "defines conventions and a set of generic components for building coupled models using the Earth System Modeling Framework (:term:`ESMF`)." 

   NWP
      Numerical Weather Prediction (NWP) takes current observations of weather and processes them with computer models to forecast the future state of the weather. 

   NWS
      The `National Weather Service <https://www.weather.gov/>`__ (NWS) is an agency of the United States government that is tasked with providing weather forecasts, warnings of hazardous weather, and other weather-related products to organizations and the public for the purposes of protection, safety, and general information. It is a part of the National Oceanic and Atmospheric Administration (NOAA) branch of the Department of Commerce.

   offline UPP
      Refers to cases where UPP is built standalone and run separately from the model.

   RAP
      `Rapid Refresh <https://rapidrefresh.noaa.gov/>`__. The continental-scale NOAA hourly-updated assimilation/modeling system operational at :term:`NCEP`. RAP covers North America and is comprised primarily of a numerical forecast model and an analysis/assimilation system to initialize that model. RAP is complemented by the higher-resolution 3km High-Resolution Rapid Refresh (:term:`HRRR`) model.

   Repository
      A central location in which files (e.g., data, code, documentation) are stored and managed. 

   RRFS
      The `Rapid Refresh Forecast System <https://gsl.noaa.gov/focus-areas/unified_forecast_system/rrfs>`__ (RRFS) is NOAA's next-generation convection-allowing, rapidly-updated, ensemble-based data assimilation and forecasting system currently scheduled for operational implementation in 2024. It is designed to run forecasts on a 3-km :term:`CONUS` domain. 

   SDF
      Suite Definition File. An external file containing information about the construction of a physics suite. It describes the schemes that are called, in which order they are called, whether they are subcycled, and whether they are assembled into groups to be called together.

   SRW
   Short-Range Weather Application
      The `Short-Range Weather Application <https://github.com/ufs-community/ufs-srweather-app>`__ is a UFS Application that targets predictions of atmospheric behavior on a limited spatial domain and on time scales from minutes out to about two days. It packages a prognostic atmospheric model (the UFS Weather Model), pre- and post-processing tools, and a community workflow.

   Spack
      `Spack <https://spack.readthedocs.io/en/latest/>`__ is a package management tool designed to support multiple versions and configurations of software on a wide variety of platforms and environments. It was designed for large supercomputing centers, where many users and application teams share common installations of software on clusters with exotic architectures. 

   spack-stack
      The `spack-stack <https://github.com/JCSDA/spack-stack>`__ is a collaborative effort between the NOAA Environmental Modeling Center (:term:`EMC`), the UCAR Joint Center for Satellite Data Assimilation (:term:`JCSDA`), and the Earth Prediction Innovation Center (:term:`EPIC`). *spack-stack* is a repository that provides a :term:`Spack`-based method for building the software stack required for numerical weather prediction (NWP) tools such as the `Unified Forecast System (UFS) <https://ufscommunity.org/>`__ and the `Joint Effort for Data assimilation Integration (JEDI) <https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/>`__ framework. *spack-stack* uses the Spack package manager along with custom Spack configuration files and Python scripts to simplify installation of the libraries required to run various applications. The *spack-stack* can be installed on a range of platforms and comes pre-configured for many systems. Users can install the necessary packages for a particular application and later add the missing packages for another application without having to rebuild the entire stack.

   UFS
      The Unified Forecast System is a community-based, coupled, comprehensive Earth modeling 
      system consisting of several applications (apps). These apps span regional to global 
      domains and sub-hourly to seasonal time scales. The UFS is designed to support the :term:`Weather Enterprise` and to be the source system for NOAA's operational numerical weather prediction applications. For more information, visit https://ufscommunity.org/.

   Updraft helicity
      Helicity measures the rotation in a storm's updraft (rising) air. Significant rotation increases the probability that the storm will produce severe weather, including tornadoes. See http://ww2010.atmos.uiuc.edu/(Gh)/guides/mtr/svr/modl/fcst/params/hel.rxml for more details on updraft helicity. 

   Weather Enterprise
      Individuals and organizations from public, private, and academic sectors that contribute to the research, development, and production of weather forecast products; primary consumers of these weather forecast products.

   Weather Model
      A prognostic model that can be used for short- and medium-range research and
      operational forecasts. It can be an atmosphere-only model or an atmospheric
      model coupled with one or more additional components, such as a wave or ocean model. The SRW App uses the `UFS Weather Model <https://github.com/ufs-community/ufs-weather-model>`__.

   Workflow
      The sequence of steps required to run an experiment from start to finish. 

   write component
      The output files written by the UFS Weather Model use an Earth System Modeling Framework (ESMF) component, referred to as the write component because the UPP cannot directly process output on the native grid types (e.g., “GFDLgrid”, “ESGgrid”). Output fields are interpolated to a write component grid before writing them to an output file. 