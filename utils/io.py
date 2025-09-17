import bioio
import xarray as xr
import scyjava as sj

def open_xarray(path: str) -> xr.DataArray:
    '''
    Use bioio (formerly aicsimageio) to open image directly as xarray
    '''
    image = bioio.BioImage(path)
    return image.xarray_data

def open_bioformats(path: str, split_channels=True) -> list:
    '''
    Use Bio-Formats Importer to 
    '''
    #quiet
    debugtools = sj.jimport('loci.common.DebugTools')
    debugtools.setRootLevel('ERROR')

    bf = sj.jimport('loci.plugins.BF')
    options = sj.jimport('loci.plugins.in.ImporterOptions')()
    options.setId(path)
    options.setOpenAllSeries(True)
    options.setSplitChannels(split_channels)
    options.setQuiet(True)
    return bf.openImagePlus(options)

def parse_channels(path: str) -> list:
    '''
    Parse channel order and colors
    '''
    #setup list
    channels = []
    
    #try parsing channel names first
    names = bioio.BioImage(path).channel_names
    for name in names:
        #green
        if 'GFP' in name:
            channels.append('G')
        #TD
        elif 'TD' in name:
            channels.append('T')
        #TODO: blue, red
    
    return channels