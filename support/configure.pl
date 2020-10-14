#  Copyright (c) 1997-2018
#  Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
#  http://www.polymake.org
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation; either version 2, or (at your option) any
#  later version: http://www.gnu.org/licenses/gpl.txt.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#-------------------------------------------------------------------------------


# This is a template for a configuration script of a bundled or standalone extension.
# Please edit it according to your needs and rename it to configure.pl
# All subroutines defined below must stay here, even if they do nothing.


# List of variables used in the ninja rules that must be modified when building this extension.
# The complete list of recognized variables is contained in the module perllib/Polymake/Configure.pm.
# Usually you will have to modify only few of them, most probably CXXFLAGS, LDFLAGS, or LIBS.
# Please put just the bare names of the variables here, without '$' prefix.

@conf_vars=qw( );


# This subroutine should augment the key sets of one or both hash maps passed by reference.
# These hash maps contain all options recognized by configure script.
# The first hash map is used for normal options like --docdir.
# The second hash map is used for trivalent options like --with-java / --without-java.
#
# The option disabling this extension completely (only applicable to bundled extensions)
# is added automatically and should not be mentioned here.
# Please specify just the bare option names, without '--' prefix.

sub allowed_options {
   my ($allowed_options, $allowed_with)=@_;
   # @$allowed_options{ qw( ACTION ) }=();
   # @$allowed_with{ qw( SOMETHING ) }=();
}


# This subroutine should print to STDERR a short explanation of specific options introduced above.

sub usage {
   # print STDERR "  --ACTION               do something extraordinary\n",
   #              "  --with-SOMETHING=PATH  installation path of SOMETHING library, if non-standard\n";
}


# This subroutine should analyze the passed options, perform whatever sanity checks,
# for example, try to compile a short test program or check a version of a required third-party library,
# and, finally, either set (some of) the variables listed in @conf_vars above, or die with a message
# explaining what was wrong and how it could be remedied.
#
# The only argument passed here is a reference to a hash map containing *all* options specified for the
# configure script, including those intended for other (bundled) extensions and the core system.
# Please be sure not to modify any options not belonging to this extension.
#
# Options specified as --without-SOMETHING are passed as (SOMETHING => ".none.") key-value pairs.
#
# You might want to use helper routines defined in the module Polymake::Configure or access the standard
# configuration variables like $InstallTop or $LDFLAGS, which are also defined in the package Polymake::Configure.
#
# Concerning make variables like CXXFLAGS or LIBS: please only assign the additional values needed
# for this extension.  They are automatically merged with other settings when this extension is being built.
# To modify make variables temporarily. e.g. in order to compile a test program, the additional values
# can be passed as optional arguments to the helper utilities `compile_test_program' and `build_test_program'.

sub proceed {
    my ($options)=@_;
}
