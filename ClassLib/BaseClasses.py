import klayout.db
from math import sqrt, cos, sin, atan2, pi, copysign
from klayout.db import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from klayout.db import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from ClassLib._PROG_SETTINGS import *

from collections import OrderedDict

class Element_Base():
    '''
    @brief: base class for simple single-layer or multi-layer elements and objects that are consisting of
            several polygons.
            metal_region polygons will be added to the design
            empty_region polygons will be erased from the background with 
            metal region polygons already added.

    '''
    def __init__(self, origin, trans_in=None, inverse=False ):
        ## MUST BE IMPLEMENTED ##
        self.connections = []       # DPoint list with possible connection points
        self.connection_edges = [] # indexes of edges that are intended to connect to other polygons
        # indexes in "self.connection_edges" where Sonnet ports
        # should be placed
        self.sonnet_port_connections = [] 
        self.angle_connections = [] #list with angle of connecting elements
        ## MUST BE IMLPEMENTED END ##
        
        self.connection_ptrs = [] # pointers to connected structures represented by their class instances
        
        self.origin = origin
        self.inverse = inverse
        self.metal_region = Region()
        self.empty_region = Region()
        self.metal_regions = {}
        self.empty_regions = {}        
        self.metal_regions["default"] = self.metal_region
        self.empty_regions["default"] = self.empty_region
        
        self.metal_region.merged_semantics = True
        self.empty_region.merged_semantics = True
        self.DCplxTrans_init = None
        self.ICplxTrans_init = None
        
        if( trans_in is not None ):
        # if( isinstance( trans_in, ICplxTrans ) ): <==== FORBIDDEN
            if( isinstance( trans_in, DCplxTrans ) ):
                self.DCplxTrans_init = trans_in
                self.ICplxTrans_init = ICplxTrans().from_dtrans( trans_in )
            elif( isinstance( trans_in, CplxTrans ) ):
                self.DCplxTrans_init = DCplxTrans().from_itrans( trans_in )
                self.ICplxTrans_init = ICplxTrans().from_trans( trans_in )
            elif( isinstance( trans_in, DTrans ) ):
                self.DCplxTrans_init = DCplxTrans( trans_in, 1 )
                self.ICplxTrans_init = ICplxTrans( Trans().from_dtrans( trans_in ), 1 )
            elif( isinstance( trans_in, Trans ) ):
                self.DCplxTrans_init = DCplxTrans( DTrans().from_itrans( trans_in ), 1 )
                self.ICplxTrans_init = ICplxTrans( trans_in, 1 )
                
        self._init_regions_trans()
    
    def init_regions( self ):
        raise NotImplementedError
    
    # first it makes trans_init displacement
    # then the rest of the trans_init
    # then displacement of the current state to the origin
    # after all, origin should be updated
    def _init_regions_trans( self ):
        self.init_regions()         # must be implemented in every subclass
        
        dr_origin = DSimplePolygon( [DPoint(0,0)] )
        if( self.DCplxTrans_init is not None ):
            # constructor trans displacement
            dCplxTrans_temp = DCplxTrans( 1,0,False, self.DCplxTrans_init.disp )
            self.make_trans( dCplxTrans_temp )
            dr_origin.transform( dCplxTrans_temp )
            
            # rest of the constructor trans functions
            dCplxTrans_temp = self.DCplxTrans_init.dup()
            dCplxTrans_temp.disp = DPoint(0,0)
            self.make_trans( dCplxTrans_temp )
            dr_origin.transform( dCplxTrans_temp )                  
            
        # translation to the old origin (self.connections are already contain proper values)
        self.make_trans( DCplxTrans( 1,0,False, self.origin ) ) # move to the origin
        self.origin += dr_origin.point( 0 )
        
    def make_trans( self, dCplxTrans ):
        if( dCplxTrans is not None ):
            iCplxTrans = ICplxTrans().from_dtrans( dCplxTrans )
            self.metal_region.transform( iCplxTrans )
            self.empty_region.transform( iCplxTrans )
            self._update_connections( dCplxTrans )
            self._update_alpha( dCplxTrans )
    
    def _update_connections( self, dCplxTrans ):       
        if( dCplxTrans is not None ):
            # the problem is, if i construct polygon with multiple points 
            # their order in poly_temp.each_point() doesn't coinside with the 
            # order of the list that was passed to the polygon constructor
            # so, when i perform transformation and try to read new values through poly_temp.each_point()
            # they values are rearranged
            # solution is: i need to create polygon for each point personally, and the initial order presists
            for i,pt in enumerate(self.connections):
                poly_temp = DSimplePolygon( [pt] )
                poly_temp.transform( dCplxTrans )
                self.connections[i] = poly_temp.point( 0 )
    
    def _update_alpha( self, dCplxTrans ):
        if( dCplxTrans is not None ):
            dCplxTrans_temp = dCplxTrans.dup()
            dCplxTrans_temp.disp = DPoint(0,0)
            
            for i,alpha in enumerate(self.angle_connections):
                poly_temp = DSimplePolygon( [DPoint( cos(alpha), sin(alpha) )] )
                poly_temp.transform( dCplxTrans_temp )
                pt = poly_temp.point( 0 )
                self.angle_connections[i] = atan2( pt.y, pt.x )
    
    def _update_origin( self, dCplxTrans ):
        if( dCplxTrans is not None ):     
            poly_temp = DSimplePolygon( [self.origin] )
            poly_temp.transform( dCplxTrans )
            self.origin = poly_temp.point( 0 )
    
    def place( self, dest, layer_i=-1, region_name=None, merge=False ):
        r_cell = None
        metal_region = None
        empty_region = None
        if( region_name == None ):
            metal_region = self.metal_region
            empty_region = self.empty_region
        else:
            if( region_name in self.metal_regions ):
                metal_region = self.metal_regions[region_name]
            else:
                metal_region = Region()
            
            if( region_name in self.empty_regions ):
                empty_region = self.empty_regions[region_name]
            else:
                empty_region = Region()
                    
        if( layer_i != -1 ): 
            r_cell = Region( dest.begin_shapes_rec( layer_i ) )        
            temp_i = dest.layout().layer( klayout.db.LayerInfo(PROGRAM.LAYER1_NUM,0) )
            r_cell += metal_region
            r_cell -= empty_region
            
            if( merge is True ):
                r_cell.merge()
                
            dest.shapes( temp_i ).insert( r_cell )
            dest.layout().clear_layer( layer_i )
            dest.layout().move_layer( temp_i, layer_i )
            dest.layout().delete_layer( temp_i )
        if( layer_i == -1 ): # dest is interpreted as instance of Region() class
            dest += metal_region
            dest -= empty_region
            if( merge is True ):
                dest.merge()
                
    def add_sonnet_port( self, connection_idx ):
        '''
        @brief: sets internal marker that during export to Sonnet
                   the port should be placed at the connection edge
                   with an index 'connection_idx' 
        '''
        print(self.connections)
        self.sonnet_port_connections.append(self.connections[connection_idx])


class Complex_Base( Element_Base ):
    def __init__( self, origin, trans_in=None ):
        super().__init__( origin,trans_in )
        self.primitives = OrderedDict() # ensures order of placing and erasing of primitives is preserved
        self._init_primitives_trans()
    
    def _init_regions_trans( self ):
        pass
    
    def make_trans( self, dCplxTrans_temp ):
        for primitive in self.primitives.values():
            primitive.make_trans( dCplxTrans_temp )
        self._update_connections( dCplxTrans_temp )
        self._update_alpha( dCplxTrans_temp )
                
    def _init_primitives_trans( self ):
        self.init_primitives()              # must be implemented in every subclass
        dr_origin = DSimplePolygon( [DPoint(0,0)] )
        if( self.DCplxTrans_init is not None ):
            # constructor trans displacement
            dCplxTrans_temp = DCplxTrans( 1,0,False, self.DCplxTrans_init.disp )
            self.make_trans( dCplxTrans_temp )
            dr_origin.transform( dCplxTrans_temp )
            
            # rest of the constructor trans functions
            dCplxTrans_temp = self.DCplxTrans_init.dup()
            dCplxTrans_temp.disp = DPoint(0,0)
            self.make_trans( dCplxTrans_temp )
            dr_origin.transform( dCplxTrans_temp )   
        
        dCplxTrans_temp = DCplxTrans( 1,0,False, self.origin )
        self.make_trans( dCplxTrans_temp ) #move to the origin
        self.origin += dr_origin.point( 0 )
        
        # FOLLOWING CYCLE GIVES WRONG INFO ABOUT FILLED AND ERASED AREAS
        for element in self.primitives.values():
            self.metal_region += element.metal_region
            self.empty_region += element.empty_region
    
    def place( self, dest, layer_i=-1, region_name = None ):
        if( layer_i != -1 ):
            r_cell = Region( dest.begin_shapes_rec( layer_i ) )
            for primitive in self.primitives.values():
                primitive.place( r_cell, region_name=region_name )
            
            temp_i = dest.layout().layer( klayout.db.LayerInfo(PROGRAM.LAYER1_NUM,0) ) 
            dest.shapes( temp_i ).insert( r_cell )
            dest.layout().clear_layer( layer_i )
            dest.layout().move_layer( temp_i, layer_i )
            dest.layout().delete_layer( temp_i )
        else:
            for primitive in self.primitives.values():
                primitive.place( dest, region_name=region_name )
            
    
    def init_primitives( self ):
        raise NotImplementedError
        
    def init_regions( self ):
        pass  
