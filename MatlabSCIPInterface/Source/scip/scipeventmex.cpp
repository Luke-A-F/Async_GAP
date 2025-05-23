/* SCIPMEX - A MATLAB MEX Interface to SCIP
 * Released Under the BSD 3-Clause License:
 *
 * Based on code by Jonathan Currie, which was based in parts on matscip.c supplied with SCIP.
 *
 * Authors:
 * - Jonathan Currie
 * - Marc Pfetsch
 */

#include <signal.h>
#include <scip/scip.h>

/* include this file last, since it redefines printf, which clashes with some definitions in SCIP */
#include "mex.h"

#ifndef HAVE_OCTAVE
/* The ut functions are private functions within Matlab; we do not need them for octave. */
#ifdef __cplusplus
/* Ctrl-C Detection */
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);
#else
extern bool utIsInterruptPending();
extern void utSetInterruptPending(bool);
#endif
#endif

/** executed when adding the event */
static
SCIP_DECL_EVENTINIT(eventInitCtrlC)
{
   /* notify SCIP to add event */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );
   return SCIP_OKAY;
}

/** executed when removing the event */
static
SCIP_DECL_EVENTEXIT(eventExitCtrlC)
{
   /* notify SCIP to drop event */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, -1) );
   return SCIP_OKAY;
}

/** executed when event occurs */
static
SCIP_DECL_EVENTEXEC(eventExecCtrlC)
{
#ifndef HAVE_OCTAVE
   /* Check for Ctrl-C */
   if (utIsInterruptPending())
   {
      utSetInterruptPending(false);
      mexPrintf("\nCtrl-C Detected. Exiting SCIP...\n\n");
      raise(SIGINT);
   }
#endif
   return SCIP_OKAY; /* always OK - otherwise we don't get intermediate answer */
}

/** add Ctrl-C event handler */
SCIP_RETCODE SCIPincludeCtrlCEventHdlr(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_EVENTHDLR* eventhdlr = NULL;

   /* create Event Handler */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, "CtrlCMatlab", "Catching Ctrl-C From Matlab", eventExecCtrlC, NULL) );

   /* Setup callbacks */
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitCtrlC) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitCtrlC) );

   return SCIP_OKAY;
}
