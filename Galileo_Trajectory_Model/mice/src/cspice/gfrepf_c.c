/*

-Procedure gfrepf_c ( GF, progress report finalization )

-Abstract

   Finish a GF progress report.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   GF
   TIME

-Keywords

   GEOMETRY
   SEARCH
   UTILITY

*/

   #include "SpiceUsr.h"
   #include "SpiceZfc.h"
   #include "SpiceZmc.h"
   #include "SpiceZst.h"

   void gfrepf_c ( void )

/*

-Brief_I/O

   VARIABLE  I/O  DESCRIPTION
   --------  ---  --------------------------------------------------
   None.

-Detailed_Input

   None.

-Detailed_Output

   None. This routine does perform console I/O when progress
   reporting is enabled.

-Parameters

   None.

-Exceptions

   1)  If an I/O error results from writing to standard output, the
       error is signaled by a routine in the call tree of this
       routine.

-Files

   None.

-Particulars

   This is one of three GF progress reporting routines that cooperate
   in order to display a report via console I/O. These routines may
   be used by SPICE-based applications as inputs to mid-level GF
   search routines.

   Developers wishing to use their own GF progress reporting routines
   must design them with the same interfaces and should assign them the
   same progress reporting roles as those of these routines.

   The GF progress reporting API routines are written to simplify
   reporting of work (such as searching for a geometric event) over a
   particular window. This is an important feature for interactive
   programs that may "go away" from the user's control for a
   considerable length of time. It allows the user to see that
   something is still going on (although maybe not too quickly).

   The three routines constituting the GF progress reporting API
   are:

      gfrepi_c  is used to prepare the reporting mechanism for a search
                pass. It is used to store the confinement window and
                progress report message prefix and suffix, and to
                initialize parameters associated with the reporting of
                the job in progress.

      gfrepu_c  is used to notify the progress reporting system that
                a specified increment of work has been completed
                since the last call to gfrepu_c or gfrepi_c, whichever
                occurred most recently.

      gfrepf_c  is used to "finish" the reporting of work (set the
                completion value to 100%.

-Examples

   The numerical results shown for this example may differ across
   platforms. The results depend on the SPICE kernels used as
   input, the compiler and supporting libraries, and the machine
   specific arithmetic implementation.

   1) This example shows how to call a mid-level GF search API that
      requires as input progress reporting routines.

      If custom progress reporting routines are available, they
      can replace gfrepi_c, gfrepu_c, and gfrepf_c in any GF API calls.

      The code example below is the first example in the header of
      gfocce_c.


      Conduct a search using default GF progress reporting
      and interrupt handling capabilities.

      The program will use console I/O to display a simple
      ASCII-based progress report.

      The program will trap keyboard interrupts (on most systems,
      generated by typing the "control C" key combination). This
      feature can be used in non-trivial applications to allow
      the application to continue after a search as been interrupted.

      The program will find occultations of the Sun by the Moon as seen
      from the center of the Earth over the month December, 2001.

      Use light time corrections to model apparent positions of Sun
      and Moon. Stellar aberration corrections are not specified
      because they don't affect occultation computations.

      We select a step size of 20 seconds, which implies we ignore
      occultation events lasting less than 20 seconds, if any exist.
      Given this step size and the length of the search interval, the
      user has time to interrupt the computation. In an interactive
      setting, the user might speed up the search by lengthening the
      step size or shortening the search interval, as long as these
      adjustments don't prevent the search from finding the correct
      solution.

      Use the meta-kernel shown below to load the required SPICE
      kernels.


         KPL/MK

         File name: gfrepf_ex1.tm

         This meta-kernel is intended to support operation of SPICE
         example programs. The kernels shown here should not be
         assumed to contain adequate or correct versions of data
         required by SPICE-based user applications.

         In order for an application to use this meta-kernel, the
         kernels referenced here must be present in the user's
         current working directory.

         The names and contents of the kernels referenced
         by this meta-kernel are as follows:

            File name                     Contents
            ---------                     --------
            de421.bsp                     Planetary ephemeris
            pck00008.tpc                  Planet orientation and
                                          radii
            naif0009.tls                  Leapseconds


         \begindata

            KERNELS_TO_LOAD = ( 'de421.bsp',
                                'pck00008.tpc',
                                'naif0009.tls'  )

         \begintext

         End of meta-kernel


      Example code begins here.


      /.
         Program gfrepf_ex1
      ./
      #include <stdio.h>
      #include "SpiceUsr.h"

      int main()
      {
         /.
         Constants
         ./
         #define  TIMFMT  "YYYY MON DD HR:MN:SC.###### ::TDB (TDB)"
         #define  CNVTOL  1.e-6
         #define  MAXWIN  200
         #define  TIMLEN  41

         /.
         Local variables
         ./
         SpiceBoolean            bail;
         SpiceBoolean            rpt;

         SpiceChar             * win0;
         SpiceChar             * win1;
         SpiceChar               begstr [ TIMLEN ];
         SpiceChar               endstr [ TIMLEN ];

         SPICEDOUBLE_CELL      ( cnfine, MAXWIN );
         SPICEDOUBLE_CELL      ( result, MAXWIN );

         SpiceDouble             et0;
         SpiceDouble             et1;
         SpiceDouble             left;
         SpiceDouble             right;

         SpiceInt                i;

         /.
         Load kernels.
         ./
         furnsh_c ( "gfrepf_ex1.tm" );

         /.
         Obtain the TDB time bounds of the confinement
         window, which is a single interval in this case.
         ./
         win0 = "2001 DEC 10 00:00:00 TDB";
         win1 = "2002 JAN 01 00:00:00 TDB";

         str2et_c ( win0, &et0 );
         str2et_c ( win1, &et1 );

         /.
         Insert the time bounds into the confinement
         window.
         ./
         wninsd_c ( et0, et1, &cnfine );

         /.
         Select a twenty-second step. We'll ignore any occultations
         lasting less than 20 seconds.
         ./
         gfsstp_c ( 20.0 );

         /.
         Turn on interrupt handling and progress reporting.
         ./
         bail = SPICETRUE;
         rpt  = SPICETRUE;

         /.
         Perform the search.
         ./
         gfocce_c ( "ANY",
                    "MOON",     "ellipsoid",  "IAU_MOON",
                    "SUN",      "ellipsoid",  "IAU_SUN",
                    "LT",       "EARTH",      CNVTOL,
                    gfstep_c,   gfrefn_c,     rpt,
                    gfrepi_c,   gfrepu_c,     gfrepf_c,
                    bail,       gfbail_c,     &cnfine,
                    &result                              );


         if ( gfbail_c() )
         {
            /.
            Clear the CSPICE interrupt indication. This is
            an essential step for programs that continue
            running after an interrupt; gfbail_c will
            continue to return SPICETRUE until this step
            has been performed.
            ./
            gfclrh_c();


            /.
            We've trapped an interrupt signal. In a realistic
            application, the program would continue operation
            from this point. In this simple example, we simply
            display a message and quit.
            ./
            printf ( "\nSearch was interrupted.\n\nThis message "
                     "was written after an interrupt signal\n"
                     "was trapped. By default, the program "
                     "would have terminated \nbefore this message "
                     "could be written.\n\n"                       );
         }
         else
         {

            if ( wncard_c(&result) == 0 )
            {
               printf ( "No occultation was found.\n" );
            }
            else
            {
               for ( i = 0;  i < wncard_c(&result);  i++ )
               {
                  /.
                  fetch and display each occultation interval.
                  ./
                  wnfetd_c ( &result, i, &left, &right );

                  timout_c ( left,  TIMFMT, TIMLEN, begstr );
                  timout_c ( right, TIMFMT, TIMLEN, endstr );

                  printf ( "Interval %d\n", (int)i );
                  printf ( "   Start time: %s\n", begstr );
                  printf ( "   Stop time:  %s\n", endstr );
               }
            }

         }

         return ( 0 );
      }


      When this program was executed on a Mac/Intel/cc/64-bit
      platform, the output was:


      Occultation/transit search 100.00% done.

      Interval 0
         Start time: 2001 DEC 14 20:10:14.195952  (TDB)
         Stop time:  2001 DEC 14 21:35:50.317994  (TDB)


      Note that the progress report has the format shown below:

         Occultation/transit search   6.02% done.

      The completion percentage was updated approximately once per
      second.

      When the program was interrupted at an arbitrary time,
      the output was:

         Occultation/transit search  13.63% done.
         Search was interrupted.

         This message was written after an interrupt signal
         was trapped. By default, the program would have terminated
         before this message could be written.

-Restrictions

   None.

-Literature_References

   None.

-Author_and_Institution

   N.J. Bachman        (JPL)
   J. Diaz del Rio     (ODC Space)
   L.S. Elson          (JPL)
   W.L. Taber          (JPL)
   I.M. Underwood      (JPL)
   E.D. Wright         (JPL)

-Version

   -CSPICE Version 1.0.1, 02-JUN-2021 (JDR)

       Edited the header to comply with NAIF standard. Added
       complete code example.

   -CSPICE Version 1.0.0, 28-FEB-2009 (NJB) (LSE) (WLT) (IMU) (EDW)

-Index_Entries

   GF finish a progress report

-&
*/

{ /* Begin gfrepf_c */

   /*
   Participate in error tracing.
   */
   if ( return_c() )
   {
      return;
   }

   chkin_c ( "gfrepf_c" );

   /*
   Let the f2c'd routine do the work.
   */
   gfrepf_ () ;


   chkout_c ( "gfrepf_c" );

} /* End gfrepf_c */
