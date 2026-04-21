#!/usr/bin/perl -w

#---------------------------------------------------------
# Utility script to summarize regression testing results
# Searches for errors in work directory
# Compares runtimes of test suites against baseline times
# Finds baseline file changes
# Copies new baseline files into place
#---------------------------------------------------------

use strict;
use Cwd qw(getcwd);
use Getopt::Long qw(GetOptions);
use File::Basename;

#------------------
# Global Variables
#------------------
my $URSA_BASELINE     = '/scratch4/NAGAPE/epic/role-epic/ursa/UPP/test_suite';
my $ORION_BASELINE    = '/work/noaa/epic/role-epic/orion/UPP';
my $HERCULES_BASELINE = '/work/noaa/epic/role-epic/hercules/UPP';
my $URSA_SUDO_PREFIX  = 'sudo su - role.epic -c';
my $MSU_SUDO_PREFIX   = 'sudo -u role-epic sh -c';
my $WORK_DIR_PATTERN  = '^work-upp-([A-Z]*)-(intel[a-z]*)$';
my $REL_LOG_DIR_PATH  = '../tests/logs';
my $RUNTIME_PATTERN   = '([a-z0-9]+)_test ([0-9:]{8}) -- baseline ([0-9:]{8})$';
my $BASELINE_PATTERN  = 'changes in results for case ([0-9a-z_]+) in (.*)$';
my $debug             = 0;
my $copy_files        = 0;

#---------------------------------------------------------
# Return the RDHPC system specific baseline directory
# Arguments: 
#   $_[0]: RDHPC system name (ursa|orion|hercules)
# Return: Path to RDHPC specific baseline directory
#---------------------------------------------------------
sub getBaselineDir()
{
  my $rv = undef;
  my $sys = shift;

  if ($sys eq "ursa")        { $rv = $URSA_BASELINE; }
  elsif ($sys eq "orion")    { $rv = $ORION_BASELINE; }
  elsif ($sys eq "hercules") { $rv = $HERCULES_BASELINE; }
  else                       { &errHand("Undefined RDHPC System: ${sys}"); }
 
  return $rv;
}

#---------------------------------------------------------
# Return the full RDHPC system specific baseline directory
# Arguments: 
#   $_[0]: RDHPC system name (ursa|orion|hercules)
#   $_[1]: Compiler
#   $_[2]: Test Suite
# Return: Path to full RDHPC specific baseline directory
#---------------------------------------------------------
sub getFullBaselineDir()
{
  my $rv =undef;
  my ($sys, $cmplr, $ts) = @_;
 
  $rv = sprintf("%s/data_out_%s/%s", &getBaselineDir($sys), $cmplr, $ts);

  &errHand("Directory DNE: ${rv} -- $!") unless(-d $rv);
  
  return $rv;
}
  

#---------------------------------------------------------
# Return the RDHPC system specific sudo command prefix
# - This will allow the script to execute as the epic
# - group user (provided you have permission). You will
# - br prompted for the password
# Arguments: 
#   $_[0]: RDHPC system name (ursa|orion|hercules)
# Return: String value of RDHPC specific sudo prefix
#---------------------------------------------------------
sub getSudoPrefix()
{
  my $rv = undef;
  my $sys = shift;
  
  if ($sys eq "ursa")                           
  { 
    $rv = $URSA_SUDO_PREFIX; 
  }
  elsif ($sys eq "hercules" or $sys eq "orion") 
  { 
    $rv = $MSU_SUDO_PREFIX; 
  }
  else                                          
  { 
    &errHand("Undefined RDHPC System: ${sys}"); 
  }
 
  return $rv;
}

#---------------------------------------------------------
# Return a search pattern to find the UPP work directory
#---------------------------------------------------------
sub getWorkDirPattern()
{
  return "${WORK_DIR_PATTERN}";
}

#---------------------------------------------------------
# Return a search pattern to find to populate runtime hash
#---------------------------------------------------------
sub getRuntimePattern()
{
  return "${RUNTIME_PATTERN}";
}

#---------------------------------------------------------
# Return a search pattern to find changed baseline files 
#---------------------------------------------------------
sub getBaselinePattern()
{
  return "${BASELINE_PATTERN}";
}

#---------------------------------------------------------
# Return the relative path to the UPP log directory
#---------------------------------------------------------
sub getRelLogDirPath
{
  return "${REL_LOG_DIR_PATH}";
}

#---------------------------------------------------------
# Change directories to provided directory
# Arguments: 
#   $_[0]: Destination directory path
#---------------------------------------------------------
sub changeDir()
{
  my $dir = shift;
  chdir($dir) or 
     &errHand("Could not change directories: $dir");
}

#---------------------------------------------------------
# Return a string with the UPP working directory
# Arguments: 
#   $_[0]: UPP Run Directory Path
# Return: String name of UPP working directory
#---------------------------------------------------------
sub getWorkDir() 
{
  my $rv = undef;
  my $dir = shift;
  my $pat = getWorkDirPattern();

  opendir(my $dh, $dir) 
    or &errHand("Could not open directory: ${dir} -- $!");
  while(readdir($dh))
  {
    $rv = $_ if(/$pat/); 
  }
  close $dh;

  return $rv;
}

#---------------------------------------------------------
# Return a string with the UPP working directory
# Arguments: 
#   $_[0]: system (lower case)
#   $_[1]: compiler
# Return: String name of log file to scrape
#---------------------------------------------------------
sub getLogFile()
{
   return "rt.log." . uc($_[0]) . "_$_[1]";
}

#---------------------------------------------------------
# Populate the runtime hash
# Arguments: 
#   $_[0]: Hash reference to runtime hash
#   $_[1]: Path to log file
#---------------------------------------------------------
sub populateRuntimes()
{
  my $rtHsh = shift;
  my $lf = shift;

  open(my $fh, '<', $lf)
    or &errHand("Unable to open the file: $lf - $!");
  while(<$fh>)
  {
    chomp;
    if ($_ =~ &getRuntimePattern())
    {
      my ($suite, $rtime, $bltime) = ($1, $2, $3);
      $rtHsh->{$suite} = [$rtime, $bltime];
    }
  }
  close($fh);
}

#---------------------------------------------------------
# Caclulate runtime deltas for each suite in seconds 
# Arguments: 
#   $_[0]: Hash reference to runtime hash
#   $_[1]: Hash reference to result hash
#---------------------------------------------------------
sub calculateRuntimes()
{
  my $rtHsh  = shift;
  my $resHsh = shift;
  foreach(keys(%{$rtHsh}))
  {
    my @run = split(':', $rtHsh->{$_}->[0]);
    my @bl  = split(':', $rtHsh->{$_}->[1]);
    my $runSecs = $run[0]*3600 + $run[1]*60 + $run[2];
    my $blSecs  = $bl[0]*3600 + $bl[1]*60 + $bl[2];
    $resHsh->{$_} = $blSecs - $runSecs;
  }
} 

#---------------------------------------------------------
# Find changed baseline files and populate baseline hash
# Arguments: 
#   $_[0]: Hash reference to baseline hash
#   $_[1]: Path to log file
#---------------------------------------------------------
sub getBaselineFiles()
{
  my $blHsh = shift;
  my $lf = shift;

  open(my $fh, '<', $lf)
    or &errHand("Unable to open the file: $lf - $!");
  while(<$fh>)
  {
    chomp;
    if ($_ =~ &getBaselinePattern())
    {
      my ($suite, $blPath) = ($1, $2);
      if ($blHsh->{$suite})
      {
        push(@{$blHsh->{$suite}}, $blPath);
      }
      else
      { 
        $blHsh->{$suite} = [$blPath]; 
      }
    }
  }
  close($fh);
}

#---------------------------------------------------------
# Create sudo command to copy baseline files
# Arguments: 
#   $_[0]: RDHPC system name (ursa|orion|hercules)
#   $_[1]: Compiler
#   $_[2]: Test Suite
#   $_[3]: Source File to be copied
# Return: String containing full sudo command to copy a
#         new baseline file
#---------------------------------------------------------
sub constructSudoCommand()
{
   my ($sys, $cmplr, $suite, $srcFile) = @_;
   my $cmd = undef;
   my $blPath = &getFullBaselineDir($sys, $cmplr, $suite);
   my $baseSrcFile = basename($srcFile);
   my $fullDestFile = sprintf("%s/%s.%s", ${blPath}, ${baseSrcFile}, uc($sys));

   # Copy new baseline file to baseline directory
   my $srcCopy = sprintf('cp %s %s/.', ${srcFile}, ${blPath});

   # Preserve a copy of the old baseline file
   my $preserveOldCopy = 'echo "New baseline file detected"';
   if (-e $fullDestFile)
   {
     $preserveOldCopy = sprintf('cp -p %s %s-old', $fullDestFile, $fullDestFile);
   }
   
   # Move the new source file into place with the full baseline file  name
   my $srcMove = sprintf('mv %s/%s %s', $blPath, $baseSrcFile, $fullDestFile);

   # Change permissions on the new source file
   my $chgPerm = sprintf('chmod 755 %s', $fullDestFile);

   # Full sudo command
   $cmd = sprintf("%s '%s;%s;%s;%s'", &getSudoPrefix($sys), $srcCopy, $preserveOldCopy, $srcMove, $chgPerm);

   return $cmd;
}

#---------------------------------------------------------
# Error Handler for script
# Arguments:
#   $_[0]: Error message
#---------------------------------------------------------
sub errHand() 
{
  printf("Error encountered: %s\n", shift); 
  exit 1;
}

#---
# Start script execution
#---
my $runDir = undef;
GetOptions(
  'run|r=s' => \$runDir,
  'debug|d' => \$debug,
  'copy|c'  => \$copy_files
) or &errHand("Usage: $0 --run <run-directory> [--copy] [--debug]");
&errHand("Error defining run directory") unless(defined($runDir));

my $workDir = &getWorkDir($runDir);
my ($rdhpc_sys, $compiler) = (lc($1), $2) if ($workDir =~ &getWorkDirPattern());
my $logFile = &getLogFile($rdhpc_sys, $compiler);

&changeDir($runDir);
chomp(my @errorArr = qx(egrep -ir "(error|fatal)" $workDir | grep -v "tee"));
&changeDir(&getRelLogDirPath);

#---
# runtime hash structure:  key -> reference array (0: run time, 1: baseline time)
# result hash structure:   key -> delta time in seconds
# baseline hash structure: key -> reference array with paths to changed baseline files
#---
my %runtimeHsh = ();
my %resultHsh = ();
my %baselineHsh = ();
&populateRuntimes(\%runtimeHsh, $logFile);
&calculateRuntimes(\%runtimeHsh, \%resultHsh);
&getBaselineFiles(\%baselineHsh, $logFile);

#---
# Print all errors in work directory
#---
print "\nChecking for errors in work directory...\n\n";
unless (@errorArr)
{
  print "#-----------------\n";
  print "There were no errors detected in the work directory\n";
  print "#-----------------\n";
}
else
{
  print "#-----------------\n";
  print "$_\n" foreach(@errorArr);
  print "#-----------------\n";
}

#---
# Print test suites timing results
#---
print "\nChecking timing results for all regresssion tests...\n\n";
print "#-----------------\n";
foreach(keys(%resultHsh))
{
  printf("%-7s is %-4d seconds %-5s time\n", 
         $_, abs($resultHsh{$_}), $resultHsh{$_} < 0 ? "OVER" : "UNDER");
}  
print "#-----------------\n";

#---
# Print test suites that have modified baseline files
# Copy modified baseline files to baseline directory
#---
print "\nChecking for baseline files that have changed...\n\n";
print "#-----------------\n";
foreach(keys(%baselineHsh))
{
  print "$_:\n";
  my $cnt = 1;
  foreach my $blf (@{$baselineHsh{$_}})
  {
    print "\t${cnt}) ${blf}\n";
    my $sudoCmd = &constructSudoCommand($rdhpc_sys, $compiler, $_, $blf);
    print "\nCommand to execute:\n---\n ${sudoCmd}\n---\n\n" if ($debug); 
    if ($copy_files) {
       if (system($sudoCmd))
       {
          &errHand("Problem encountered while copying baseline file: ${sudoCmd}")
       }
       print "\t   ***New baseline file copied successfully***\n";
    }
    $cnt++;
  }
} 
print "No baseline file changes\n" unless(keys(%baselineHsh));
print "#-----------------\n";

#---------------------------------------------------------
# Strictly for debugging values stored in the variables
#---------------------------------------------------------
if ($debug) 
{
  print "Run Directory   ==> ${runDir}\n";
  print "Work Directory  ==> ${workDir}\n";
  print "RDHPC System    ==> ${rdhpc_sys}\n";
  print "System Compiler ==> ${compiler}\n";
  print "Log File        ==> ${logFile}\n";
  foreach(@errorArr)
  {
    print "$_\n";
  }
  foreach(keys(%runtimeHsh))
  {
    print "$_:\n";
    print "---Run time:      $runtimeHsh{$_}->[0]\n";
    print "---Baseline time: $runtimeHsh{$_}->[1]\n";
  }
  foreach(keys(%resultHsh))
  {
    print "$_:\n";
    print "---Delta: $resultHsh{$_}\n";
  }
  foreach(keys(%baselineHsh))
  {
    print "$_:\n";
    foreach my $blf (@{$baselineHsh{$_}})
    {
      print "==>${blf}\n";
    }
  }
}
