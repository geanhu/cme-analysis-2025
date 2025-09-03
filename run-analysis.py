import argparse
import sys

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
    args = parser.parse_args()

    #Launch GUI or enter workflow directly
    if args.gui:
        gui()
    elif args.workflow is None:
        print('Must set --workflow if -gui is disabled')
        sys.exit(1)
    else:
        workflow()


def gui():
    '''
    TODO: Launch GUI
    '''
    return

def workflow():
    '''
    TODO: Launch workflow
    '''
    return

if __name__ == "__main__":
    main()