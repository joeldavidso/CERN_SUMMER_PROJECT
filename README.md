These files are split into two sections,
      - The event display code which consists of:
      	- EVENTDISPLAYS (where the event displays are saved)
	- ATLASSTYLE (if you need it)
	- project_EDv2.C 
	  - This is the main event display code, it is all explained in the code and there should only be 1 section you change to use it.
	- project_EDvP.C
	  - This is the test/proof of concept of the polar event display view, it works exactly the same as the main event display in use.

      - The Processing code:
      	- PLOTS
	  - This is where the drawn hists are saved
      	- project_Hv2.h
	  - This is a small header file that contains some variable that may be shared between the processing code and the drawing code (although its easy enough to not need this)
	- project_Pv3.C 
	  - This is the main processing code which then saves all the histograms to a file to be later read and drawn by the drawing code
	- project_Mv3.C 
	  - This is the main drawing code which reads the processed hists and draws them, it has options to draw normalise hists and to sleect which ones to draw individually