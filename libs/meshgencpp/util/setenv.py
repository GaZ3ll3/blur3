# ============================================================================
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import os
import getopt

def write_environment_variable(csh_handle,bash_handle,var,value):
    csh_handle.write( 'setenv %s "%s"\n' % (var.upper(),value))
    bash_handle.write('export %s="%s"\n' % (var.upper(),value))

def write_shortcut_alias(csh_handle,bash_handle,var,value):
    csh_handle.write( 'alias %s "%s"\n' % (var,value))
    bash_handle.write('alias %s="%s"\n' % (var,value))


# Paths
if os.path.exists("../qhull"):
    qhull_path = os.path.abspath("../qhull")
else:
    print ""
    qhull_tmp = raw_input(' Enter path to qhull package: ')
    if os.path.exists(qhull_tmp):
        qhull_path = os.path.abspath(qhull_tmp)
        print ""
        print " Found qhull directory: ",qhull_path
        print ""
    else:
        print ""
        print " Could not find qhull directory: ",qhull_tmp
        print ""
        print " Make sure that qhull has been properly installed"
        print "   and that you provided the proper path."
        print ""
        print " Once you are certain qhull is installed and you"
        print "   know the proper path, please re-run this script by typing:"
        print ""
        print " $ python util/setenv.py"
        print ""
        exit();
    
meshgen_path = os.path.abspath(os.curdir)
outfile_base="setenv"
    
# Open output files
csh_file = open(os.path.join(meshgen_path,".".join((outfile_base,"csh"))),'w')
bash_file = open(os.path.join(meshgen_path,".".join((outfile_base,"bash"))),'w')

# Write out boiler plate
boiler_plate = ("# MeshGenC++ environment settings\n")
csh_file.write(boiler_plate)
bash_file.write(boiler_plate)

# Write out variables
python_path = "${PYTHONPATH}"
matlab_path = "${MATLABPATH}"

meshgenP = "".join((meshgen_path,"/viz/python"))
meshgenM = "".join((meshgen_path,"/viz/matlab"))

print meshgenP
print meshgenM

python_path = ":".join((meshgenP,python_path))
matlab_path = ":".join((meshgenM,matlab_path))

print ""
print "The following variables will be set:"
print "  MESHGENCPP = %s" % meshgen_path
print "  QHULL = %s" % qhull_path
write_environment_variable(csh_file,bash_file,"MESHGENCPP",meshgen_path)
write_environment_variable(csh_file,bash_file,"QHULL",qhull_path)

print "  PYTHONPATH = %s" % python_path
write_environment_variable(csh_file,bash_file,"PYTHONPATH",python_path)

print "  MATLABPATH = %s" % matlab_path
write_environment_variable(csh_file,bash_file,"MATLABPATH",matlab_path)

pmesh2_command = "python $MESHGENCPP/viz/python/plotmesh2.py"
print "  plotmesh2 = %s" %pmesh2_command
write_shortcut_alias(csh_file,bash_file,"plotmesh2",pmesh2_command)
print ""


# Close output files
csh_file.close()
bash_file.close()
