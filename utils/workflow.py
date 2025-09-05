class workflow:
    def __init__(self):
        self.steps = []
        self.files = []

class file_filter:
    def __init__(self, 
                 allowed_ext = ['.tif', '.tiff', '.nd2', '.czi'], 
                 banned_names = ['DIC'],
                 banned_start = ['.']):
        self.allowed = allowed_ext
        self.banned = banned_names
        self.banned_start = banned_start
    
    def filter(self, files: list):
        filtered = []
        for file in files:
            filter_out = False
            #Filter out files with banned phrases in name
            for name in self.banned:
                if name in file:
                    filter_out = True
            #Filter out files starting with banned start characters
            for name in self.banned_start:
                if file[0] == name:
                    filter_out = True
            if filter_out:
                continue

            #Filter for only allowed extensions
            filter_in = False
            for ext in self.allowed:
                if file.endswith(ext):
                    filter_in = True
            if filter_in:
                filtered.append(file)
            else:
                continue
        return filtered

    

def launch_workflow(workflow_dict: dict, directory: str, ij):
    new_workflow = workflow()

    #TODO: parse steps

    #Filter list of files to process
    if workflow_dict['custom_filter']:
        directory_file_filter = file_filter(
            workflow_dict['allowed_extensions'],
            workflow_dict['banned_names'],
            workflow_dict['banned_start']
        )
    else:
        directory_file_filter = file_filter()
    #TODO: Find all files in directory -> run through file filter

    return 0

def main():
    return

if __name__ == "main":
    main()