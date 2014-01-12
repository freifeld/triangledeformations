#!/usr/bin/env python
"""
Created on Sat Jun 29 16:55:38 2013

Author: Oren Freifeld
Email: freifeld@dam.brown.edu
"""

try:
    try: # for ipython version < 0.13
        from IPython.Shell import IPShellEmbed
        ipshell = IPShellEmbed()  
    except ImportError:  # for 0.13 <= ipython version < 1.0.0
        try:     
            from IPython.frontend.terminal.embed import InteractiveShellEmbed            
            ipshell = InteractiveShellEmbed( banner1 = 'Dropping into IPython',
                        exit_msg = 'Leaving IPython, back to program.')
                          
        except TypeError: # for 1.0.0 <= ipython version 
            from IPython.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed( banner1 = 'Dropping into IPython',
                        exit_msg = 'Leaving IPython, back to program.')
    
except:
    print """
    Failed importing ipython embedded shell. Pitty. 
    This is useful, but not a must. 
    """        
    ipshell=None   
