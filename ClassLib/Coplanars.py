import klayout.db
from math import sqrt, cos, sin, atan2, pi, copysign, tan
from klayout.db import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from klayout.db import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from ClassLib.BaseClasses import *

class CPWParameters:
    '''
    Base CPW parameters to inherit from

    Parameters:
        width: float
            Width of the internal conductor
        gap: float
            Gap between the internal conductor and the ground planes
    '''
  def __init__(self, width, gap):
    self.width = width
    self.gap = gap
    self.b = 2*gap + width

class CPW( Element_Base ):
    '''
    Base class representing a single coplanar waveguide

    Parameters:
        width: float
            Width of the internal conductor. If None, it will use CPWParameters.width
        gap: float
            Gap between the internal conductor and the ground planes. If None, it will use CPWParameters.width
        start: DPoint
            Center aligned point, determines the start point of the coplanar segment
        stop: DPoint
            Center aligned point, determines the end point of the coplanar segment
        gndWidth: float
            Width of ground plane to be drawn
        trans_in:
            Klayout transformation to be made
        cpw_params: CPWParameters 
            Base CPWParameters you can inherit from. If None, it will use width and gap.
    '''
    def __init__(self, width=None, gap=None, start=DPoint(0,0), end=DPoint(0,0), gndWidth=-1, trans_in=None, cpw_params=None ):
        if( cpw_params  is None ):
            self.width = width
            self.gap = gap
            self.b = 2*gap + width
        else:
            self.width = cpw_params.width
            self.gap = cpw_params.gap
            self.b = 2*self.gap + self.width
            
        self.gndWidth = gndWidth
        self.end = end
        self.start = start
        self.dr = end - start
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_regions( self ):
        self.connections = [DPoint(0,0),self.dr]
        self.start = DPoint(0,0)
        self.end = self.start + self.dr
        alpha = atan2( self.dr.y, self.dr.x )
        self.angle_connections = [alpha,alpha]
        alpha_trans = ICplxTrans().from_dtrans( DCplxTrans( 1,alpha*180/pi,False, self.start ) )
        
        metal_poly = DSimplePolygon( [DPoint(0,-self.width/2),
                                                       DPoint(self.dr.abs(),-self.width/2), 
                                                       DPoint(self.dr.abs(),self.width/2),
                                                       DPoint(0,self.width/2)] )
        self.connection_edges = [3,1]
        self.metal_region.insert( klayout.db.SimplePolygon().from_dpoly( metal_poly ) )
        if( self.gap != 0 ):
            self.empty_region.insert( klayout.db.Box( Point().from_dpoint(DPoint(0,self.width/2)), Point().from_dpoint(DPoint( self.dr.abs(), self.width/2 + self.gap )) ) )
            self.empty_region.insert( klayout.db.Box( Point().from_dpoint(DPoint(0,-self.width/2-self.gap)), Point().from_dpoint(DPoint( self.dr.abs(), -self.width/2 )) ) )
        self.metal_region.transform( alpha_trans )
        self.empty_region.transform( alpha_trans )
        
        
class Air_bridges( Element_Base ):
    '''
    Base class representing airbridges

    Parameters:
        N_air_bridges: float
            Number of airbridges to be drawn
        width: float
            Width of the internal conductor. If None, it will use CPWParameters.width
        gap: float
            Gap between the internal conductor and the ground planes. If None, it will use CPWParameters.width
        start: DPoint
            Center aligned point, determines the start point of the coplanar segment
        stop: DPoint
            Center aligned point, determines the end point of the coplanar segment
        gndWidth: float
            Width of ground plane to be drawn
        trans_in:
            Klayout transformation to be made
        cpw_params: CPWParameters 
            Base CPWParameters you can inherit from. If None, it will use width and gap.
    '''
    def __init__(self,  N_air_bridges, width=None, gap=None, start=DPoint(0,0), end=DPoint(0,0), gndWidth=-1, trans_in=None, cpw_params=None ):
        if( cpw_params  is None ):
            self.width = width
            self.gap = gap
            self.b = 2*gap + width
        else:
            self.width = cpw_params.width
            self.gap = cpw_params.gap
            self.b = 2*self.gap + self.width
            
        self.gndWidth = gndWidth
        self.N_air_bridges = N_air_bridges
        self.end = end
        self.start = start
        self.dr = end - start
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_regions( self ):
        self.connections = [DPoint(0,0),self.dr]
        self.start = DPoint(0,0)
        self.end = self.start + self.dr
        self.L0 = self.start.distance(self.end) / (self.N_air_bridges +1)
        alpha = atan2( self.dr.y, self.dr.x )
        self.angle_connections = [alpha,alpha]
        alpha_trans = ICplxTrans().from_dtrans( DCplxTrans( 1,alpha*180/pi,False, self.start ) )
        for i in range(self.N_air_bridges):
            self.air_bridge = klayout.db.DBox( DPoint((i+1)*self.L0, -self.b/2 - 2e3), DPoint((i+1)*self.L0+2e3, self.b/2 + 2e3) ) 
            self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge ) )
            
            self.air_bridge_2 = klayout.db.DBox( DPoint((i+1)*self.L0-self.b/4+1e3, -self.b/2 - 2e3), DPoint((i+1)*self.L0+1e3+self.b/4, -self.b/2 - 4e3) ) 
            self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_2 ) )
            
            self.air_bridge_3 = klayout.db.DBox( DPoint((i+1)*self.L0-self.b/4+1e3, self.b/2 + 2e3), DPoint((i+1)*self.L0+1e3+self.b/4, self.b/2 + 4e3) ) 
            self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_3 ) )
       
        self.metal_region.transform( alpha_trans )


class CPW_arc( Element_Base ):
    '''
    Base class representing a single coplanar waveguide arc

    Parameters:
        Z0: CPW or CPWParameters
            Parameters to inherit width and gap from
        start: DPoint
            Center aligned point, determines the start point of the coplanar segment
        R: float
            Radius of the arc
        delta_alpha: float
            Angle of the arc
        gndWidth: float
            Width of the ground to be drawn
        trans_in: Transformation
            KLayout transformation to be processed during execution
    '''
    def __init__(self, Z0, start, R, delta_alpha, gndWidth = -1, trans_in=None ):
        self.R = R
        self.start = start
        self.center = start + DPoint( 0,self.R )
        self.end = self.center + DPoint( sin(delta_alpha), -cos(delta_alpha) )*self.R
        self.dr = self.end - self.start
        
        self.width = Z0.width
        self.gap = Z0.gap
        
        self.delta_alpha = delta_alpha
        self.alpha_start = 0
        self.alpha_end = self.delta_alpha

        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.center = self.connections[2]
        
        #print("End coords:", self.dr, self.end)
        
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def _get_solid_arc( self, center, R, width, alpha_start, alpha_end, n_inner, n_outer ):
        pts = []
#        print(alpha_start/pi, alpha_end/pi, cos( alpha_start ), cos( alpha_end ),
#                         sin(alpha_start), sin(alpha_end))
                         
        if alpha_end > alpha_start:
            alpha_start = alpha_start - 1e-3 
            alpha_end = alpha_end + 1e-3
        else:
            alpha_start = alpha_start + 1e-3
            alpha_end = alpha_end - 1e-3
        
        d_alpha_inner = (alpha_end - alpha_start)/(n_inner - 1)
        d_alpha_outer = -d_alpha_inner
        
        
#        print("Center:", center)
        
        for i in range(0, n_inner):
            alpha = alpha_start + d_alpha_inner*i
            pts.append(center + DPoint(cos(alpha), sin(alpha))*(R - width/2))
        for i in range(0, n_outer):
            alpha = alpha_end + d_alpha_outer*i
            pts.append(center + DPoint(cos(alpha), sin(alpha))*(R + width/2))
#        print("Points:", pts[:n_inner],"\n       ", pts[n_inner:], "\n")
        return DSimplePolygon(pts)
        
    def init_regions(self):
        self.connections = [DPoint(0, 0), self.dr, DPoint(0, self.R)]
        self.angle_connections = [self.alpha_start, self.alpha_end]
        self.start = DPoint(0, 0)
        self.end = self.dr
        self.center = DPoint(0, self.R)
        
        n_inner = 200
        n_outer = 200
        
        metal_arc = self._get_solid_arc(self.center, self.R, self.width, 
                    self.alpha_start - pi/2, self.alpha_end - pi/2, n_inner, n_outer)  
        self.connection_edges = [n_inner+n_outer,n_inner]
        
        empty_arc1 = self._get_solid_arc(self.center, self.R - (self.width + self.gap)/2, 
                    self.gap, self.alpha_start - pi/2, self.alpha_end - pi/2, n_inner, n_outer)  
        
        empty_arc2 = self._get_solid_arc(self.center, self.R + (self.width + self.gap)/2, 
                    self.gap, self.alpha_start - pi/2, self.alpha_end - pi/2, n_inner, n_outer)  
        
        self.metal_region.insert(SimplePolygon().from_dpoly(metal_arc))
        self.empty_region.insert(SimplePolygon().from_dpoly(empty_arc1))
        self.empty_region.insert(SimplePolygon().from_dpoly(empty_arc2))
        
        
class CPW_arc_2( Element_Base ):
# Just a small change compared to CPW_arc, now alpha-start is an input parameter so that 
# we can start with a non-horizontal segment 
    def __init__(self, Z0, start, R, delta_alpha, alpha_start, gndWidth = -1, trans_in=None ):
        self.R = R
        self.start = start
        self.delta_alpha = delta_alpha
        self.alpha_start = alpha_start
        self.alpha_end = self.delta_alpha + self.alpha_start
        
        
        self.center = self.start - self.R * DPoint( cos(self.alpha_start), sin(self.alpha_start) )
        self.end = self.center + DPoint( cos(self.alpha_end), sin(self.alpha_end) )*self.R
        self.dr = self.end - self.start
        
        self.width = Z0.width
        self.gap = Z0.gap

        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.center = self.connections[2]
        
        #print("End coords:", self.dr, self.end)
        
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def _get_solid_arc( self, center, R, width, alpha_start, alpha_end, n_inner, n_outer ):
        pts = []
#        print(alpha_start/pi, alpha_end/pi, cos( alpha_start ), cos( alpha_end ),
#                         sin(alpha_start), sin(alpha_end))
                         
        if alpha_end > alpha_start:
            alpha_start = alpha_start - 1e-3 
            alpha_end = alpha_end + 1e-3
        else:
            alpha_start = alpha_start + 1e-3
            alpha_end = alpha_end - 1e-3
        
        d_alpha_inner = (alpha_end - alpha_start)/(n_inner - 1)
        d_alpha_outer = -d_alpha_inner
        
        
#        print("Center:", center)
        
        for i in range(0, n_inner):
            alpha = alpha_start + d_alpha_inner*i
            pts.append(center + DPoint(cos(alpha), sin(alpha))*(R - width/2))
        for i in range(0, n_outer):
            alpha = alpha_end + d_alpha_outer*i
            pts.append(center + DPoint(cos(alpha), sin(alpha))*(R + width/2))
#        print("Points:", pts[:n_inner],"\n       ", pts[n_inner:], "\n")
        return DSimplePolygon(pts)
        
    def init_regions(self):
        self.connections = [DPoint(0, 0), self.dr,- self.R * DPoint( cos(self.alpha_start), sin(self.alpha_start) ) ]
        self.angle_connections = [self.alpha_start, self.alpha_end]
        self.start = DPoint(0, 0)
        self.end = self.dr
        self.center =  self.start - self.R * DPoint( cos(self.alpha_start), sin(self.alpha_start) )
        
        n_inner = 200
        n_outer = 200
        
        metal_arc = self._get_solid_arc(self.center, self.R, self.width, 
                    self.alpha_start , self.alpha_end , n_inner, n_outer)  
        self.connection_edges = [n_inner+n_outer,n_inner]
        
        empty_arc1 = self._get_solid_arc(self.center, self.R - (self.width + self.gap)/2, 
                    self.gap, self.alpha_start, self.alpha_end, n_inner, n_outer)  
        
        empty_arc2 = self._get_solid_arc(self.center, self.R + (self.width + self.gap)/2, 
                    self.gap, self.alpha_start, self.alpha_end, n_inner, n_outer)  
        
        self.metal_region.insert(SimplePolygon().from_dpoly(metal_arc))
        self.empty_region.insert(SimplePolygon().from_dpoly(empty_arc1))
        self.empty_region.insert(SimplePolygon().from_dpoly(empty_arc2))
        
        
class Air_bridges_arc( Element_Base ):

    '''
    Base class representing airbridges along a curved element

    Parameters:
    '''
    def __init__(self, Z0, start, R, delta_alpha, gndWidth = -1, trans_in=None ):
        self.R = R
        self.start = start
        self.center = start + DPoint( 0,self.R )
        self.end = self.center + DPoint( sin(delta_alpha), -cos(delta_alpha) )*self.R
        self.dr = self.end - self.start
        
        self.width = Z0.width
        self.gap = Z0.gap
        self.b = 2*self.gap + self.width
        
        self.delta_alpha = delta_alpha
        self.alpha_start = 0
        self.alpha_end = self.delta_alpha

        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.center = self.connections[2]
        
        #print("End coords:", self.dr, self.end)
        
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

        
    def init_regions(self):
        self.connections = [DPoint(0, 0), self.dr, DPoint(0, self.R)]
        self.angle_connections = [self.alpha_start, self.alpha_end]
        self.start = DPoint(0, 0)
        self.end = self.dr
        self.center = DPoint(0, self.R)
        
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
        Delta=[0,0,0]
        Vec=[0,0,0]
        VecOrt=[0,0,0]
        for i, alpha in enumerate([self.alpha_start-pi/2, self.alpha_start+self.delta_alpha/2-pi/2, self.alpha_end-pi/2]):
            Delta[i]= (self.b/2 + 2e3)*DPoint(cos(alpha), sin(alpha)) + 1e3* DPoint(-sin(alpha), cos(alpha))
            Vec[i]= DPoint(cos(alpha), sin(alpha))
            VecOrt[i] = DPoint(-sin(alpha), cos(alpha))
            
        self.air_bridge_10 = klayout.db.DBox(self.start + Delta[0], self.start - Delta[0]) 
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_10 ) )
        self.air_bridge_11 = klayout.db.DBox(self.start + Delta[0] + (self.b/4-1e3)*VecOrt[0], self.start + Delta[0] + (-self.b/4-1e3)*VecOrt[0]+ 2e3*Vec[0])  
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_11 )) 
        self.air_bridge_12 = klayout.db.DBox(self.start - Delta[0] - (self.b/4-1e3)*VecOrt[0], self.start - Delta[0] + (self.b/4+1e3)*VecOrt[0]- 2e3*Vec[0])  
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_12 ) )
        
        self.medium = self.center + DPoint( sin(self.delta_alpha/2), -cos(self.delta_alpha/2) )*self.R
       
        self.air_bridge_20 = klayout.db.DBox(self.medium + Delta[1], self.medium -Delta[1] ) 
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_20 ) )
        self.air_bridge_21 = klayout.db.DBox(self.medium + Delta[1] + (self.b/4-1e3)*VecOrt[1], self.medium + Delta[1] + (-self.b/4-1e3)*VecOrt[1]+ 2e3*Vec[1])  
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_21 )) 
        self.air_bridge_22 = klayout.db.DBox(self.medium - Delta[1] - (self.b/4-1e3)*VecOrt[1], self.medium - Delta[1] + (self.b/4+1e3)*VecOrt[1]- 2e3*Vec[1])  
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_22 ) )
        
        self.air_bridge_30 = klayout.db.DBox(self.end + Delta[2], self.end - Delta[2] ) 
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_30 ) )
        self.air_bridge_31 = klayout.db.DBox(self.end + Delta[2] + (self.b/4-1e3)*VecOrt[2], self.end + Delta[2] + (-self.b/4-1e3)*VecOrt[2]+ 2e3*Vec[2])  
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_31 )) 
        self.air_bridge_32 = klayout.db.DBox(self.end - Delta[2] - (self.b/4-1e3)*VecOrt[2], self.end - Delta[2] + (self.b/4+1e3)*VecOrt[2]- 2e3*Vec[2])  
        self.metal_region.insert( klayout.db.Box().from_dbox( self.air_bridge_32 ) )
        

class CPW2CPW( Element_Base ):
    '''
    Base class representing a smooth connection between two coplanar waveguides.

    Parameters:
        Z0: CPW or CPWParameters
            First CPW parameters
        Z1: CPW or CPWParameters
            Second CPW parameters
        start: DPoint
            Center aligned point, determines the start point of the coplanar segment
        end: Dpoint
            Center aligned point, determines the end point of the coplanar segment
        trans_in: Transformation
            Klayout transformation to be executed
    '''
    def __init__(self, Z0, Z1, start, end, trans_in = None):
        self.Z0 = Z0
        self.Z1 = Z1
        self.start = start
        self.end = end
        self.dr = self.end - self.start
        super().__init__(start, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_regions(self):
        self.connections = [DPoint(0, 0), DPoint(self.dr.abs(), 0)]
        self.angle_connections = [0, 0]
        alpha = atan2(self.dr.y, self.dr.x)
        self.angle_connections = [alpha,alpha]
        alpha_trans = DCplxTrans(1, alpha*180/pi, False, 0, 0)
        
        m_poly = DSimplePolygon([DPoint(0,-self.Z0.width/2), DPoint(self.dr.abs(), -self.Z1.width/2), 
                                    DPoint(self.dr.abs(),self.Z1.width/2), DPoint(0,self.Z0.width/2)] )
        e_poly1 = DSimplePolygon([DPoint(0,-self.Z0.b/2), DPoint(self.dr.abs(), -self.Z1.b/2), 
                                    DPoint(self.dr.abs(),-self.Z1.width/2), DPoint(0,-self.Z0.width/2)] )
        e_poly2 = DSimplePolygon([DPoint(0,self.Z0.b/2), DPoint(self.dr.abs(),self.Z1.b/2), 
                                    DPoint(self.dr.abs(),self.Z1.width/2), DPoint(0,self.Z0.width/2)] )
        
        m_poly.transform(alpha_trans)
        e_poly1.transform(alpha_trans)
        e_poly2.transform(alpha_trans)
        
        self.metal_region.insert(SimplePolygon.from_dpoly(m_poly))
        self.empty_region.insert(SimplePolygon.from_dpoly(e_poly1))
        self.empty_region.insert(SimplePolygon.from_dpoly(e_poly2))
  
class Path_RSRS( Complex_Base ):
    def __init__( self, Z0, start, L1, L2, R1, R2, trans_in=None ):
        self.Z0 = Z0
        self.L1 = L1
        self.L2 = L2
        self.R1 = R1
        self.R2 = R2
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
        
    def init_primitives( self ):
        self.arc1 = CPW_arc( self.Z0, DPoint(0,0), self.R1, pi/4 )
        self.cop1 = CPW( self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end + DPoint( sin(pi/4), sin(pi/4) )*self.L1 )
        self.arc2 = CPW_arc( self.Z0, self.cop1.end, self.R2, pi/4, trans_in=DCplxTrans( 1,45,False,0,0 ) )
        self.cop2 = CPW( self.Z0.width, self.Z0.gap,  self.arc2.end, self.arc2.end + DPoint( 0,self.L2 ) )
        
        self.connections = [self.arc1.start,self.cop2.end]
        self.angle_connections = [self.arc1.alpha_start, self.cop2.alpha_end]
        
        self.primitives = {"arc1":self.arc1, "cop1":self.cop1, "arc2":self.arc2, "cop2":self.cop2}
        
class Path_RS(Complex_Base):
    def __init__(self, Z0, start, end, trans_in=None):
        self.Z0 = Z0
        self.end = end
        self.dr = end - start
        super().__init__(start, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]

    def init_primitives(self):
        L = abs(abs(self.dr.y) - abs(self.dr.x))
        r = min(abs(self.dr.y), abs(self.dr.x))
        if(abs(self.dr.y) > abs(self.dr.x)):
            self.arc1 = CPW_arc(self.Z0, DPoint(0, 0), copysign(
                r, self.dr.y), copysign(pi/2, self.dr.x*self.dr.y))
            self.cop1 = CPW(self.Z0.width, self.Z0.gap, self.arc1.end,
                            self.arc1.end + DPoint(0, copysign(L, self.dr.y)))
            self.connections = [self.arc1.start, self.cop1.end]
            self.angle_connections = [self.arc1.alpha_start, self.cop1.alpha_end]
        else:
            self.cop1 = CPW(self.Z0.width, self.Z0.gap, DPoint(
                0, 0), DPoint(0, 0) + DPoint(copysign(L, self.dr.x), 0))
            self.arc1 = CPW_arc(self.Z0, self.cop1.end, copysign(
                r, self.dr.y), copysign(pi/2, self.dr.y*self.dr.x))
            self.connections = [ self.cop1.start, self.arc1.end]
            self.angle_connections = [ self.cop1.alpha_start,self.arc1.alpha_end]
            
        self.primitives = {"arc1": self.arc1, "cop1": self.cop1}

        
class Coil_type_1( Complex_Base ):
# On this version, the initial CPW section is horizontal, 
# in order to have a horizontal initial CPW we create Coil_type_2
    
    def __init__( self, Z0, start, L1, r, L2, trans_in=None ):
        self.Z0 = Z0
        self.L1 = L1
        self.r = r
        self.L2 = L2
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_primitives( self ):
            self.cop1 = CPW( self.Z0.width, self.Z0.gap, DPoint(0,0), DPoint( self.L1,0 ) )
            self.arc1 = CPW_arc( self.Z0, self.cop1.end, -self.r, -pi )
            self.cop2 = CPW( self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end - DPoint( self.L2,0 ) )
            self.arc2 = CPW_arc( self.Z0, self.cop2.end, -self.r, pi )
            
            self.connections = [self.cop1.start,self.arc2.end]
            self.angle_connections = [self.cop1.alpha_start,self.arc2.alpha_end]
            self.primitives = {"cop1":self.cop1,"arc1":self.arc1,"cop2":self.cop2,"arc2":self.arc2}
            
            
class Coil_type_2( Complex_Base ):
# On this version, the initial CPW section is vertical, 

    
    def __init__( self, Z0, start, L1, r, L2, initial_direction = "down", trans_in=None ):
        self.Z0 = Z0
        self.L1 = L1
        self.r = r
        self.L2 = L2
        self.initial_direction = initial_direction
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
        
    def init_primitives( self ):
        
        if self.initial_direction == "down":   
            self.cop1 = CPW( self.Z0.width, self.Z0.gap, DPoint(0,0), DPoint(0, -self.L1 ) )
            self.arc1 = CPW_arc_2( self.Z0, self.cop1.end, -self.r, pi, 0 )
            self.cop2 = CPW( self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end + DPoint( 0, self.L2 ) )
            self.arc2 = CPW_arc_2( self.Z0, self.cop2.end, -self.r, -pi, 0 )
        
        else:
            self.cop1 = CPW( self.Z0.width, self.Z0.gap, DPoint(0,0), DPoint(0, +self.L1 ) )
            self.arc1 = CPW_arc_2( self.Z0, self.cop1.end, -self.r, -pi, 0 )
            self.cop2 = CPW( self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end + DPoint( 0,- self.L2 ) )
            self.arc2 = CPW_arc_2( self.Z0, self.cop2.end, -self.r, pi, 0 )
            
        self.connections = [self.cop1.start,self.arc2.end]
        self.angle_connections = [self.cop1.alpha_start,self.arc2.alpha_end]
        self.primitives = {"cop1":self.cop1,"arc1":self.arc1,"cop2":self.cop2,"arc2":self.arc2}
    

    
class Coil_air_bridges( Complex_Base ):
    def __init__( self, Z0, start, L1, r, L2, N_air_bridges, trans_in=None ):
        self.Z0 = Z0
        self.L1 = L1
        self.r = r
        self.L2 = L2
        self.N_air_bridges = N_air_bridges
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_primitives( self ):
            self.cop1 = Air_bridges(self.N_air_bridges, self.Z0.width, self.Z0.gap, DPoint(0,0), DPoint( self.L1,0 ) )
            self.arc1 = Air_bridges_arc( self.Z0, self.cop1.end, -self.r, -pi )
            self.cop2 = Air_bridges(self.N_air_bridges, self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end - DPoint( self.L2,0 ))
            self.arc2 = Air_bridges_arc( self.Z0, self.cop2.end, -self.r, pi )
            
            self.connections = [self.cop1.start,self.arc2.end]
            self.angle_connections = [self.cop1.alpha_start,self.arc2.alpha_end]
            self.primitives = {"cop1":self.cop1,"arc1":self.arc1,"cop2":self.cop2,"arc2":self.arc2}
            
            
            
            
class Coil_air_bridges_type_2( Complex_Base ):
# Contraru to Coil_air_bridges, this one only places the air bridges in the straight lines, we place nothing on 
# the turns of the coil.
    def __init__( self, Z0, start, L1, r, L2, N_air_bridges, trans_in=None ):
        self.Z0 = Z0
        self.L1 = L1
        self.r = r
        self.L2 = L2
        self.N_air_bridges = N_air_bridges
        super().__init__( start, trans_in )
        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_primitives( self ):
            self.cop1 = Air_bridges(self.N_air_bridges, self.Z0.width, self.Z0.gap, DPoint(0,0), DPoint( self.L1,0 ) )
            self.cop2 = Air_bridges(self.N_air_bridges, self.Z0.width, self.Z0.gap, self.arc1.end, self.arc1.end - DPoint( self.L2,0 ))
            
            self.connections = [self.cop1.start,self.arc2.end]
            self.angle_connections = [self.cop1.alpha_start,self.arc2.alpha_end]
            self.primitives = {"cop1":self.cop1,"arc1":self.arc1,"cop2":self.cop2,"arc2":self.arc2}
            
            
            
            
from collections import Counter

class CPW_RL_Path(Complex_Base):

    def __init__(self, origin, shape, cpw_parameters, turn_radiuses, 
        segment_lengths, turn_angles, trans_in = None):
        '''
        A piecewise-linear coplanar waveguide with rounded turns.
        
        Segment lengths are treated as the lengths of the segments of
        a line with turn_radiuses = 0. Changing turning radiuses
        will not alter the position of the end of the line.
        
        TODO: 180 deg turns
        
        Parameters:
            origin: DPoint
                The point where the line should start 
            shape: string
                String in format "RLLRL" where an R means a turn 
                and an L means a straight part
            cpw_parameters: CPWParameters or array-like
                Parameters of the CPW or an array-like with parameters 
                for each peace (R or L)
            turn_radiuses: float or array-like
                Radius of the turns or an array-like with radiuses for 
                each turn 
            segment_lengths: array-like
                Lengths of the straight parts of the equivalent 
                piecewise-linear line with no corner rounding
            turn_angles: array-like
                Angles for each turn of the line
                !!! 180 turns are not yet supported !!!
            trans_in: DTrans
                Transformation of the line as a whole            
        '''
    
        self._shape_string = shape
        self._N_elements = len(shape)
        self._shape_string_counter = Counter(shape)
        
        self._N_turns = self._shape_string_counter['R']
        
        if hasattr(cpw_parameters, "__len__"):
            if len(cpw_parameters) != self._N_elements:
                raise ValueError("CPW parameters dimension mismatch")
            self._cpw_parameters = cpw_parameters
        else:
            self._cpw_parameters = [cpw_parameters]*self._N_elements
        
        if hasattr(turn_radiuses, "__len__"):
            if len(turn_radiuses) !=  self._N_turns:
                raise ValueError("Turn raduises dimension mismatch")
            self._turn_radiuses = turn_radiuses
        else:
            self._turn_radiuses = [turn_radiuses]* self._N_turns
    
        self._segment_lengths = segment_lengths
        
        if hasattr(turn_angles, "__len__"):
            if len(turn_angles) !=  self._N_turns:
                raise ValueError("Turn angles dimension mismatch")
            self._turn_angles = turn_angles
        else:
            self._turn_angles = [turn_angles]* self._N_turns
        
        super().__init__(origin, trans_in)
        self.start = self.connections[0]
        self.end = self.connections[1]
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
    
        
    def init_primitives(self):

        R_index = 0
        L_index = 0
        origin = DPoint(0,0)
        
        prev_primitive_end = origin
        prev_primitive_end_angle = 0
        
        for i, symbol in enumerate(self._shape_string):
                              
            if( symbol == 'R' ):
                if( self._turn_angles[R_index] > 0 ):
                    turn_radius = self._turn_radiuses[R_index]
                else:
                    turn_radius =-self._turn_radiuses[R_index]
                                        
                cpw_arc = CPW_arc(self._cpw_parameters[i], prev_primitive_end, 
                          turn_radius, self._turn_angles[R_index], 
                            trans_in = DCplxTrans(1, prev_primitive_end_angle*180/pi, False, 0, 0))
                
                self.primitives["arc_"+str(R_index)] = cpw_arc
                R_index += 1    
                
            elif symbol == 'L':
            
                # Turns are reducing segments' lengths so as if there were no roundings at all
                
                # next 'R' segment if exists
                if( i+1 < self._N_elements
                    and self._shape_string[i+1] == 'R' 
                    and abs(self._turn_angles[R_index]) < pi ):
                            coeff = abs(tan(self._turn_angles[R_index]/2))
                            self._segment_lengths[L_index] -= self._turn_radiuses[R_index]*coeff
                # previous 'R' segment if exists
                if( i-1 > 0
                    and self._shape_string[i-1] == 'R' 
                    and abs(self._turn_angles[R_index-1]) < pi ):    
                        coeff = abs(tan(self._turn_angles[R_index-1]/2))
                        self._segment_lengths[L_index] -= self._turn_radiuses[R_index-1]*coeff
                        
                #if( self._segment_lengths[L_index] < 0 ):
                #    print(self._segment_lengths[L_index])
                #    print("CPW_RL_Path warning: segment length is less than zero")
                #    print("L_index = {}".format(L_index))

                    
                cpw = CPW(self._cpw_parameters[i].width, self._cpw_parameters[i].gap,
                        prev_primitive_end, prev_primitive_end + DPoint(self._segment_lengths[L_index], 0),
                            trans_in=DCplxTrans(1, prev_primitive_end_angle*180/pi, False, 0, 0))
                            
                self.primitives["cpw_"+str(L_index)] = cpw
                L_index += 1 
            
            primitive = list(self.primitives.values())[i]
            prev_primitive_end = primitive.end
            prev_primitive_end_angle = primitive.alpha_end
                
        self.connections = [origin, list(self.primitives.values())[-1].end]
        self.angle_connections = [0, list(self.primitives.values())[-1].alpha_end]
        
    def get_total_length(self):
        return sum(self._segment_lengths) + \
            sum([abs(R*alpha) for R, alpha in zip(self._turn_radiuses, self._turn_angles)])


            
