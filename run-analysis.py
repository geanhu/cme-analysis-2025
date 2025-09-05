import argparse
import sys
import imagej
import json

from gui.gui_utils import launch_gui
from utils.workflow import launch_workflow

def main():
    #Parse arguments
    parser = argparse.ArgumentParser(
        description="Launch analysis workflow"
    )
    parser.add_argument(
        '-gui',
        type=bool,
        default=True,
        help="Disable option to skip launching GUI (-w and -d must be set if GUI disabled)"
    )
    parser.add_argument(
        '-w', '--workflow',
        type=str,
        default=None,
        help="Path to .json file storing past workflow settings. Create and save workflows through the GUI first."
    )
    parser.add_argument(
        '-d', '--directory',
        type=str,
        default=None,
        help="Path to directory with image files to process."
    )
    parser.add_argument(
        '-ij', '--imagej',
        type=str,
        default=None,
        help="Path to local installation of ImageJ (or Fiji) to launch processes from."
    )
    parser.add_argument(
        '-v', '--version',
        type=str,
        default='2.16.0',
        help='Version of Fiji to launch if no local installation is supplied. Overwritten by -ij option.'
    )
    args = parser.parse_args()

    #Launch GUI to fetch settings
    if args.gui:
        launch_gui()
    elif args.workflow is None:
        print('Must set --workflow if -gui is disabled')
        sys.exit(1)
    
    #Initialize ImageJ
    if args.imagej is None:
        ij = imagej.init(f'sc.fiji:fiji:{args.version}')
    else:
        ij = imagej.init(args.imagej)

    #Launch workflow
    with open(args.workflow, 'r') as file:
        workflow = json.load(file)
    exit_code = launch_workflow(workflow, args.directory, ij)

    #Conclude
    if exit_code == 0:
        print("Finished workflow. Exiting session.")
        sys.exit(0)

if __name__ == "__main__":
    main()