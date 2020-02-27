import klayout.db
from klayout.db import Region
from ClassLib import PROGRAM

from collections import OrderedDict

class Chip_Design:
    """ @brief:     inherit this class for working on a chip design
                    and override draw() method where other drawing
                    methods should be called from
                    call show() to draw everything
        @params:    str cell_name - name of a cell, e.g. 'testScript'
    """
    def __init__(self, cell_name):
        # getting main references of the application
        app = klayout.db.Application.instance()
        mw = app.main_window()
        self.lv = mw.current_view()
        self.cv = None
        self.cell = None
        self.region_ph = Region()
        self.region_el = Region()
        
        #this insures that lv and cv are valid objects
        if(self.lv == None):
            self.cv = mw.create_layout(1)
            self.lv = mw.current_view()
        else:
            self.cv = self.lv.active_cellview()

        # find or create the desired by programmer cell and layer
        layout = self.cv.layout()
        layout.dbu = 0.001
        if(layout.has_cell(cell_name)):
            self.cell = layout.cell(cell_name)
        else:
            self.cell = layout.create_cell(cell_name)
        
        info = klayout.db.LayerInfo(1, 0)
        info2 = klayout.db.LayerInfo(2, 0)
        self.layer_ph = layout.layer(info)  # photoresist layer
        self.layer_el = layout.layer(info2)  # e-beam lithography layer

        # clear this cell and layer
        self.cell.clear()

        # setting layout view  
        self.lv.select_cell(self.cell.cell_index(), 0)
        self.lv.add_missing_layers()

        # design parameters that were passed to the last
        # self.draw(...) call are stored here as ordered dict
        self.design_pars = OrderedDict()
    
    # Call other methods drawing parts of the design from here
    def draw(self, design_params=None):
        '''
        @brief: Purely virtual base-class method that is ought to be
                implemented in child classes.
                Responsible for calling functions that draw separate
                objects.

                Must be started with self.deisgn_pars = design_params
        @params:
                design_params - OrderedDict()
                    dictionary that contains design parameters and
                    used by other drawing routines
        '''
        raise NotImplementedError

    # Call this m
    def show(self, design_params=None):
        self.__transfer_reg2cell()

    def __transfer_reg2cell(self):
        # this too methods assumes that all previous drawing
        # functions are placing their object on regions
        # in order to avoid extensive copying of the polygons
        # to/from cell.shapes during the logic operations on
        # polygons
        self.cell.shapes(self.layer_ph).insert(self.region_ph)
        self.cell.shapes(self.layer_el).insert(self.region_el)
        self.lv.zoom_fit()
    
    # Erases everything outside the box
    def crop(self, box, layer=None):
        if layer is None:
            self.__erase_in_layer(self.layer_ph, box)
            self.__erase_in_layer(self.layer_el, box)
        else:
            self.__erase_in_layer(layer, box)

    # Erases everything outside the box in a layer
    def __erase_in_layer(self, layer, box):
        reg_l = self._reg_from_layer(layer)
        box_reg = Region(box)
        reg_l &= box_reg

        # r_cell = r_cell - inverse_region

    def _reg_from_layer(self, layer):
        if layer == self.layer_el:
            return self.region_el
        elif layer == self.layer_ph:
            return self.region_ph
        else:
            return None

    # Save your design as GDS-II
    def save_as_gds2(self, filename):
        slo = klayout.db.SaveLayoutOptions()
        slo.format = 'GDS2'
        slo.gds2_libname = 'LIB'
        slo.gds2_max_cellname_length = 32000
        slo.gds2_max_vertex_count = 8000
        slo.gds2_write_timestamps = True
        slo.select_all_layers()
        self.lv.save_as(self.cell.cell_index(), filename, slo)