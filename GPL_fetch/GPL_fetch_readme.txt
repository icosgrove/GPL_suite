GPL_Fetch Process 1: Manual-GPL Detection Pairing Input Guide

Once your excel file containing logged detections is entered into GPL_fetch.m script under 'user input' and your parameters are set in GPL_fetch_parameter_input.m a plot window will show up allowing you yo view the manual detection in the GPL window. For user entry:

- Any number >0 will record that number GPL detection as a single call
	- if the number is > number of calls on the window, redo
- Selecting 0 will open for special cases, this should be selected for ANY case beyond pairing a single call and moving on. (Including pairing a single call in combination with other things).

Special Cases (case sensitive)
- r: Regular entry
- m: Multiple GPL detection span single detection
- e: GPL miss
- a: Adhoc detection (save other GPL detections in the window)
- s: Switch manual detection to another (better) call
- f: False log entry
- rp1: Reprocess: GPL detection missed part of the call
- rp2: Reprocess: GPL window spans call, but the call is weak
- fn: Separate call in window is good, but GPL missed it: Flag window
- new: Other call type in window: Flag window
- c: Terminate special case entry and move on to the next window
- x: User flags an error that they made, in the case that the ability to go back and fix it isn't implemented
- k: Skip: Acknowledge that the GPL detector did flag the call, but it is poor and not worth saving, and the 	rest of the window is also poor and not worth saving
- und: Call cannot be determined due to window crossing
- b: Breakout: If you enter into a special case mode and wish to break out of it
- end: End session (All pairings,adhocs are automatically saved to output file after each entry). Only enter this AFTER you have completed the current pairing. (Ex: Pair Call 7 in special case mode, when done enter 'end'. The next session will begin on Call 8)
- set: Enter settings mode to adjust brightness and contrast
	- when adjusting the B/C for the GPL contour window, I recommend increasing the contrast up around 	 	1000-3000, then adding brightness as needed (in increments of 10, 10-60+)