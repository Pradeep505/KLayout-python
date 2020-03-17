import klayout.db
from math import sqrt, cos, sin, atan2, pi, copysign
from klayout.db import Point,DPoint,DSimplePolygon,SimplePolygon, DPolygon, Polygon,  Region
from klayout.db import Trans, DTrans, CplxTrans, DCplxTrans, ICplxTrans

from importlib import reload

from .BaseClasses import Complex_Base
from .Coplanars import CPW_RL_Path, CPW
from .Coplanars import Coil_type_1, Coil_type_2, Coil_air_bridges, CPW_arc, CPW_arc_2 # backward compatibility TODO: delete classes
from .Capacitors import Connector_large_capacitor


class CPWResonator():
    
    _c = 299792458.0
    
    def __init__(self, origin, cpw_parameters, turn_radius, frequency, ɛ, 
        wavelength_fraction = 1/4, coupling_length = 200e3, offset_length = 200e3,
        meander_periods = 4, neck_length = 100e3, end_primitive = None, trans_in = None):
        '''
        A CPW resonator of wavelength fraction with automatic length calculation
        from given frequency.
        
        It's also possible to attach a primitive to the end of the resonator which
        should provide a get_phase_shift(   ) method to calculate it's effective length
        as if it was just a straight CPW piece.
        
        Parameters:
            origin: DPoint
                The point where the resonator's coupling tail ends
                and the meander starts
            cpw_parameters: CPWParameters
                CPW parameters for the resonator
            turn_radius: float
                Raduis of the resonator's meander turns 
            frequency: float, GHz
                The frequency that the resonator will have
            ɛ: float
                Substrate relative permittivity
            wavelength_fraction: float
                The fraction of the wavelength the resonator's 
                fundamental mode has at resonance. This parameter
                is not checked for consistensy with the resonator's 
                boundaries
            coupling_length: float
                The length of the segment parallel to the feedline
            offset_length: float
                The distance away from the feedline before meandering
            meander_periods: int
                The number of periods in the meandering part
            neck_length: float
                The length of the straight part before the uncoupled end
                of the resonator
            end_primitive: Element_Base object
                Element that will be attached to the end of the resonator
        '''
        self._origin = origin
        self._cpw_parameters = cpw_parameters
        self._turn_radius = turn_radius
        self._coupling_length = coupling_length
        self._offset_length = offset_length
        self._meander_periods = meander_periods
        self._neck_length = neck_length
        self._frequency = frequency
        self._wavelength_fraction = wavelength_fraction
        self._ɛ = ɛ
        self._end_primitive = end_primitive
        self._length = self._calculate_total_length()
        self._trans_in = trans_in
        
        
    def place(self, cell, layer):
        N, R = self._meander_periods, self._turn_radius
        
        meander_length = (self._length - self._coupling_length \
            - self._offset_length - 2*2*pi*R/4 \
            -N*2*pi*R - 2*pi*R/2 + R - 2*pi*R/4 - self._neck_length)/(2*N+3/2)
        shape = "LRLRL"+N*"RLRL"+"RLRL"
        segment_lengths = [self._coupling_length + R, self._offset_length + 2*R]\
             + [meander_length+R] + 2*N*[meander_length] \
             + [meander_length/2, self._neck_length+R]
        turn_angles = [pi/2, pi/2] + N*[-pi,pi]+[-pi, pi/2]
        
        line = CPW_RL_Path(self._origin, shape, 
            self._cpw_parameters, self._turn_radius, 
                segment_lengths, turn_angles, self._trans_in)
        line.place(cell, layer)
        
        self.start = line.connections[0]
        self.end = line.connections[1]
        self.alpha_start = line.angle_connections[0]
        self.alpha_end = line.angle_connections[1]
        
        
    def _calculate_total_length(self):
        length = self._c/self._frequency/(sqrt(self._ɛ/2+0.5))/1e9*self._wavelength_fraction
        # if self._end_primitive is not None:
        #     claw_phase_shift = self._claw.get_phase_shift(self._frequency)
        #     claw_effective_length = 1/180*claw_phase_shift/self._frequency/1e9*self._c/2/sqrt(self._ɛ/2+0.5)/2
        #     length -= claw_effective_length
        return length*1e9 # in nm

class CPWResonator2(Complex_Base):
    # neck doesn't stand out and instead spreads along the resonator
    _c = 299792458.0
    
    def __init__(self, origin, cpw_parameters, turn_radius, frequency, ɛ, 
        wavelength_fraction = 1/4, coupling_length = 200e3,
        meander_periods = 4, neck_length = 100e3, no_neck=False, extra_neck_length=0, end_primitive=None, trans_in = None):
        """
        A CPW resonator of wavelength fraction with automatic length calculation
        from given frequency.

        Modification: with a neck along the resonator

        It is also possible to attach a primitive to the end of the resonator which
        should provide a get_phase_shift(   ) method to calculate its effective length
        as if it was just a straight CPW piece.

        Parameters:
            origin: DPoint
                The point where the resonator's coupling tail starts
            cpw_parameters: CPWParameters
                CPW parameters for the resonator
            turn_radius: float
                Raduis of the resonator's meander turns
            frequency: float, GHz
                The frequency that the resonator will have
            ɛ: float
                Substrate relative permittivity
            wavelength_fraction: float
                The fraction of the wavelength the resonator's
                fundamental mode has at resonance. This parameter
                is not checked for consistency with the resonator's
                boundaries
            coupling_length: float
                The length of the segment parallel to the feedline
            meander_periods: int
                The number of periods in the meandering part
            neck_length: float
                The length of the straight part before the uncoupled end
                of the resonator
            no_neck: bool
                If no_neck is true then the neck is absent
            extra_neck_length: float
                The length of the segment perpendicular to the neck
        """
        self.Z = cpw_parameters
        self.R = turn_radius
        self.cpl_L = coupling_length
        self.N = meander_periods
        self.neck = neck_length
        self.freq = frequency
        self.wavelength_fraction = wavelength_fraction
        self.eps = ɛ
        self.no_neck = no_neck
        self.extra_neck_length = extra_neck_length
        self.len = self._calculate_total_length()

        self.line = None # RL path that represents meandr
        self.open_end = end_primitive # open-ended coplanar at the and of the resonator

        super().__init__(origin, trans_in)

        self.start = self.line.connections[0]
        self.end = self.line.connections[1]
        self.alpha_start = self.line.angle_connections[0]
        self.alpha_end = self.line.angle_connections[1]
        
    def init_primitives(self):
        origin = DPoint(0, 0)
        total_neck_length =  pi * self.R / 2 + self.neck if not self.no_neck else 0
        meander_length = (self.len - self.cpl_L - 2 * (1 + self.N) * pi * self.R -
                         3 * self.R - total_neck_length - self.extra_neck_length) / (2 * self.N + 3 / 2)
        shape = "LRL" + self.N * "RLRL" + "RL" + ("RL" if not self.no_neck else "")
        segment_lengths = [self.cpl_L, meander_length] + 2 * self.N * [meander_length] + 
                        [meander_length + 3 * self.R + self.extra_neck_length] + ([self.neck] if not self.no_neck else [])
        turn_angles = [pi] + self.N * [-pi, pi]+[-pi] + ([-pi/2] if not self.no_neck else [])
        
        self.line = CPW_RL_Path(origin, shape, self.Z, self.R, segment_lengths, turn_angles)
        self.primitives['line'] = self.line

        if self.open_end == None:
            self.open_end = CPW(0, self.Z.b/2, self.line.end, self.line.end +
                                 (DPoint(0, -self.Z.b) if not self.no_neck else DPoint(-self.Z.b, 0)))
        self.primitives['open_end'] = self.open_end
        
    def _calculate_total_length(self):
        length = self._c / self.freq / sqrt(self.eps/2 + 0.5) / 1e9 * self.wavelength_fraction
        return length * 1e9 # in nm

class EMResonator_TL2Qbit_worm( Complex_Base ):
    '''
        A CPW Resonator with three sections. The three sections are a straight line, a coupler, and a worm with an open end.
        
        Parameters:
            Z0: CPW
                A CPW object defining the width and gap of the resonator
            start: DPoint
                The starting point of the resonator
            L_coupling: float
                The length of the segment parallel to the feedline
            L_1: float
                Length of each one of the straight sections of the worm
            r: float
                Radius of the resonator's meander turns
            L_2: float
                Length of the open end segment
            N: int
                The number of periods in the meandering part
            trans_in: 
                Klayout transformation to be made, such as displacement, rotation etc
        '''
    def __init__( self, Z0, start, L_coupling, L1, r, L2, N, trans_in=None ):
        self.Z0 = Z0
        self.L_coupling = L_coupling
        self.L1 = L1
        self.r = r
        self.L2 = L2

        self.N = N
        self.primitives_gnd = {}
        super().__init__( start, trans_in )

        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_primitives( self ):
        name = "coil0"
        setattr(self, name, Coil_type_1(self.Z0, DPoint(0,0), self.L_coupling, self.r, self.L1) )
        self.primitives[name] = getattr( self, name )      
        # coils filling        
        for i in range(self.N):
            name = "coil" + str(i+1)
            setattr( self, name, Coil_type_1( self.Z0, DPoint( -self.L1 + self.L_coupling, -(i+1)*(4*self.r) ), self.L1, self.r, self.L1 ) )
            self.primitives[name] = getattr( self, name )
            
        # draw the "tail"
        self.arc_tail = CPW_arc( self.Z0, self.primitives["coil" + str(self.N)].end, -self.L1/2, -pi/2 )
        self.cop_tail = CPW( self.Z0.width, self.Z0.gap, self.arc_tail.end, self.arc_tail.end - DPoint( 0,self.L2 ) )
        self.cop_open_end = CPW( 0, self.Z0.b/2, self.cop_tail.end, self.cop_tail.end - DPoint(0,self.Z0.b) )
        self.primitives["arc_tail"] = self.arc_tail
        self.primitives["cop_tail"] = self.cop_tail
        self.primitives["cop_open_end"] = self.cop_open_end
        
        self.connections = [DPoint(0,0), self.cop_tail.end]
        self.angle_connections = [0,self.cop_tail.alpha_end]
        
class Air_bridges_TL2Qbit_worm( EMResonator_TL2Qbit_worm ):
    '''
        A CPW Resonator which inherits from EMResonator_TL2Qbit_worm.
        This class also draws airbridges along the length of the resonator.
        
        Parameters:
            Z0: CPW
                A CPW object defining the width and gap of the resonator
            start: DPoint
                The starting point of the resonator
            L_coupling: float
                The length of the segment parallel to the feedline
            L_1: float
                Length of each one of the straight sections of the worm
            r: float
                Radius of the resonator's meander turns
            L_2: float
                Length of the open end segment
            N: int
                The number of periods in the meandering part
            N_air_bridges: int
                Number of airbridges to setup along the length of the resonator
            trans_in: 
                Klayout transformation to be made, such as displacement, rotation etc
        '''
    def __init__( self, Z0, start, L_coupling, L1, r, L2, N, N_air_bridges, trans_in=None ):
        self.N = N
        self.N_air_bridges = N_air_bridges
        super().__init__(Z0, start, L_coupling, L1, r, L2, N)

        self.L0 = self.L1 / (self.N_air_bridges +1)
        
    def init_primitives( self ):
        name = "coil0"
        setattr(self, name, Coil_air_bridges(self.Z0, DPoint(0,0), self.L_coupling, self.r, self.L1, self.N_air_bridges) )
        self.primitives[name] = getattr( self, name )      
        # coils filling        
        for i in range(self.N):
            name = "coil" + str(i+1)
            setattr( self, name, Coil_air_bridges( self.Z0, DPoint( -self.L1 + self.L_coupling, 
                    -(i+1)*(4*self.r) ), self.L1, self.r, self.L1, self.N_air_bridges) )
            self.primitives[name] = getattr( self, name )
               
        self.connections = [DPoint(0,0), self.primitives["coil" + str(self.N)].end]
        self.angle_connections = [0,self.primitives["coil" + str(self.N)].alpha_end]

class EMResonator_TL2Qbit_worm2( EMResonator_TL2Qbit_worm ):
    '''
        A CPW Resonator with three sections. The three sections are a straight line, a coupler, and a worm with an open end.        
        Parameters:
            Z0: CPW
                A CPW object defining the width and gap of the resonator
            start: DPoint
                The starting point of the resonator
            L_coupling: float
                The length of the segment parallel to the feedline
            L_1: float
                Length of each one of the straight sections of the worm
            r: float
                Radius of the resonator's meander turns
            L_2: float
                Length of the open end segment
            N: int
                The number of periods in the meandering part
            trans_in: 
                Klayout transformation to be made, such as displacement, rotation etc
        '''
    def __init__( self, Z0, start, L_coupling, L1, r, L2, N, trans_in=None ):
        super()__init__( Z0, start, L_coupling, L1, r, L2, N, trans_in )
        
    def init_primitives( self ):
        self.arc1 = CPW_arc( self.Z0, DPoint(0,0), -self.r, pi/2 )
        self.primitives["arc1"] = self.arc1
        
        # making coil
        name = "coil0"
        setattr(self,name, Coil_type_1(self.Z0, DPoint(0,0), self.L_coupling, self.r, self.L1) )
        self.primitives[name] = getattr( self, name )      
        # coils filling        
        for i in range(self.N):
            name = "coil" + str(i+1)
            setattr( self, name, Coil_type_1( self.Z0, DPoint( -self.L1 + self.L_coupling,
                                             -(i+1)*(4*self.r) ), self.L1, self.r, self.L1 ) )
            self.primitives[name] = getattr( self, name )
            
        # draw the "tail"
        self.arc_tail = CPW_arc( self.Z0, self.primitives["coil" + str(self.N)].end, -self.L1/2, -pi/2 )
        self.cop_tail = CPW( self.Z0.width, self.Z0.gap, self.arc_tail.end, self.arc_tail.end - DPoint( 0,self.L2 ) )
        self.cop_open_end = CPW( 0, self.Z0.b/2, self.cop_tail.end, self.cop_tail.end - DPoint(0,self.Z0.b) )
        self.primitives["arc_tail"] = self.arc_tail
        self.primitives["cop_tail"] = self.cop_tail
        self.primitives["cop_open_end"] = self.cop_open_end
        
        self.connections = [DPoint(0,0), self.cop_tail.end]
        self.angle_connections = [0,self.cop_tail.alpha_end]



class Resonator_Test( EMResonator_TL2Qbit_worm ):
    def __init__( self, Z0, start, L_coupling, L1, r, L2, N, angle, trans_in=None ):
        '''
        A CPW Resonator with three sections. The three sections are a straight line, a coupler, and a worm with an open end.
        This version can change the angle of the open end.
        
        Parameters:
            Z0: CPW
                A CPW object defining the width and gap of the resonator
            start: DPoint
                The starting point of the resonator
            L_coupling: float
                The length of the segment parallel to the feedline
            r: float
                Radius of the resonator's meander turns
            L_2: float
                Length of the open end segment
            N: int
                The number of periods in the meandering part
            N_air_bridges: int
                Number of airbridges to setup along the length of the resonator
            trans_in: 
                Klayout transformation to be made, such as displacement, rotation etc
        '''
        super().__init__( Z0, start, L_coupling, L1, r, L2, N, trans_in )
        self.angle = angle
        
    def init_primitives( self ):
        name = "coil0"
        setattr(self, name, Coil_type_1(self.Z0, DPoint(0,0), self.L_coupling, self.r, self.L1) )
        self.primitives[name] = getattr( self, name )      
        # coils filling        
        for i in range(self.N):
            name = "coil" + str(i+1)
            setattr( self, name, Coil_type_1( self.Z0, DPoint( -self.L1 + self.L_coupling, -(i+1)*(4*self.r) ),
                                             self.L1, self.r, self.L1 ) )
            self.primitives[name] = getattr( self, name )
            
        # draw the "tail"
        self.arc_tail = CPW_arc( self.Z0, self.primitives["coil" + str(self.N)].end, -self.L1/8, self.angle )
        self.cop_tail = CPW( self.Z0.width, self.Z0.gap, self.arc_tail.end, self.arc_tail.end - DPoint( 0,self.L2 ) )
        self.cop_open_end = CPW( 0, self.Z0.b/2, self.cop_tail.end, self.cop_tail.end - DPoint(0,self.Z0.b) )
        self.primitives["arc_tail"] = self.arc_tail
        self.primitives["cop_tail"] = self.cop_tail
        self.primitives["cop_open_end"] = self.cop_open_end
        
        self.connections = [DPoint(0,0), self.cop_tail.end]
        self.angle_connections = [0,self.cop_tail.alpha_end]
        
class Three_Part_Resonator( EMResonator_TL2Qbit_worm ):

    '''
        A CPW Resonator with three sections. The three sections are a straight line, a coupler, and a worm with an open end.
        This version can change the angle of the open end.
        
        Parameters:
            Z0: CPW
                A CPW object defining the width and gap of the resonator
            start: DPoint
                The starting point of the resonator
            L_coupling: float
                The length of the segment parallel to the feedline
            r: float
                Radius of the resonator's meander turns
            L_2: float
                Length of the open end segment
            N: int
                The number of periods in the meandering part
            angle: float
                The angle of the open end
            trans_in: 
                Klayout transformation to be made, such as displacement, rotation etc
        '''
    def __init__( self, Z0, start, L_tail, L_coupling, L1, r, L2, N, angle, trans_in=None ):
        self.L_tail = L_tail
        self.angle = angle
        super().__init__( Z0, start, L_coupling, L1, r, L2, N, trans_in )

    def init_primitives( self ):
        name = "coil0"
        setattr(self, name, Coil_type_1(self.Z0, DPoint(0,0), self.L_coupling, self.r, self.L1) )
        self.primitives[name] = getattr( self, name )      
        # coils filling        
        for i in range(self.N):
            name = "coil" + str(i+1)
            setattr( self, name, Coil_type_1(self.Z0, DPoint( -self.L1 + self.L_coupling, -(i+1)*(4*self.r) ),
                                            self.L1, self.r, self.L1 ) )
            self.primitives[name] = getattr( self, name )
            
        # draw the "tail"
        self.arc_tail = CPW_arc( self.Z0, self.primitives["coil" + str(self.N)].end, -self.L1/8, self.angle )
        self.cop_tail = CPW( self.Z0.width, self.Z0.gap, self.arc_tail.end, self.arc_tail.end - DPoint( 0,self.L2 ) )
        self.cop_open_end = CPW( 0, self.Z0.b/2, self.cop_tail.end, self.cop_tail.end - DPoint(0,self.Z0.b) )
        self.shorted_end_arc = CPW_arc(self.Z0, DPoint(0, 0), -self.r, pi/2)
        self.shorted_end = CPW( self.Z0.width, self.Z0.gap, self.shorted_end_arc.end,
                                self.shorted_end_arc.end + DPoint( 0, -self.L_tail ) )
        self.primitives["arc_tail"] = self.arc_tail
        self.primitives["cop_tail"] = self.cop_tail
        self.primitives["cop_open_end"] = self.cop_open_end
        self.primitives["shorted_end_arc"] = self.shorted_end_arc
        self.primitives["shorted_end"] = self.shorted_end
        
        self.connections = [DPoint(0,0), self.cop_tail.end]
        self.angle_connections = [0,self.cop_tail.alpha_end, self.shorted_end_arc.alpha_end]

class EMResonator_worm( EMResonator_TL2Qbit_worm ):
    '''
        A CPW Resonator with three sections. The three sections are a straight line, a coupler, and a worm with an open end.

        Parameters:
            Z0: CPW
                A CPW object defining the width and gap of the resonator
            start: DPoint
                The starting point of the resonator
            L_coupling: float
                The length of the segment parallel to the feedline
            r: float
                Radius of the resonator's meander turns
            L_2: float
                Length of the open end segment
            N: int
                The number of periods in the meandering part
            trans_in: 
                Klayout transformation to be made, such as displacement, rotation etc
        '''
    def __init__( self, Z0, start, L_coupling, L1, r, L2, N, trans_in=None ):
        super().__init__(Z0, start, L_coupling, L1, r, L2, N, trans_in)
        
    def init_primitives( self ):
        name = "coil0"
        setattr(self, name, Coil_type_1(self.Z0, DPoint(0,0), self.L_coupling, self.r, self.L1) )
        self.primitives[name] = getattr( self, name )      
        # coils filling        
        for i in range(self.N):
            name = "coil" + str(i+1)
            setattr( self, name, Coil_type_1( self.Z0, DPoint( -self.L1 + self.L_coupling, -(i+1)*(4*self.r) ),
                                            self.L1, self.r, self.L1 ) )
            self.primitives[name] = getattr( self, name )
            
        # draw the "tail"
        self.arc_tail = CPW_arc( self.Z0, self.primitives["coil" + str(self.N)].end, -self.L1/2, -pi/2 )
        self.cop_tail = CPW( self.Z0.width, self.Z0.gap, self.arc_tail.end, self.arc_tail.end - DPoint( 0,self.L2 ) )
        self.cop_open_end = CPW( 0, self.Z0.b/2, DPoint(0, 0), DPoint(0, 0) - DPoint(self.Z0.b, 0) )
        self.primitives["arc_tail"] = self.arc_tail
        self.primitives["cop_tail"] = self.cop_tail
        self.primitives["cop_open_end"] = self.cop_open_end
        
        self.connections = [DPoint(0,0), self.cop_tail.end]
        self.angle_connections = [0,self.cop_tail.alpha_end]
        
        
class Purcell_filtered_resonator( Complex_Base ):
    def __init__( self, Z0, start, L_sep, L_connector, L_qubit_coupling, Sep_connector, L_wigling,
                 r, L2_r, L2_f, N_wigling, gap_2, spacing, d_1, d_2, N, trans_in = None ):
        self.Z0 = Z0
        # Width of the capacitor
        self.d_1 = d_1
        # Space between 2 pads of capacitor 
        self.d_2 = d_2
        # Number of connections between the 2 pads 
        self.N = N
        # Width of the line connecting the capacitor to circuit
        self.width = self.Z0.width
        # Gap near connecting wire
        self.gap_1 = self.Z0.gap
        # Gap around the capacitor 
        self.gap_2 = gap_2 
        # Spacint between lines of capacitor 
        self.spacing = spacing
        # Distance between transmission line and capacitor
        self.L_sep = L_sep
        # Lenght of the wigling, common to filter and readout cavity
        self.L_wigling = L_wigling
        # Aditional length after wigling to further tune resonator frequency
        self.L2_r = L2_r
        self.r = r
        # Additional length after wigling to further tune Purcell filter frequency
        self.L2_f = L2_f
        # Number of wigling,common to filter and readout cavity
        self.N_wigling = N_wigling
        # Distance between the readout cavity and the Purcell filter
        self.L_connector = L_connector
        # Separation between the 2 sides of the connection between Purcell filter and readout
        self.Sep_connector = Sep_connector
        # Extra length of the resonator to which we couple the qubit
        self.L_qubit_coupling = L_qubit_coupling
        

        self.primitives_gnd = {}
        super().__init__( start, trans_in )

        self.start = self.connections[0]
        self.end = self.connections[-1]
        self.dr = self.end - self.start
        self.alpha_start = self.angle_connections[0]
        self.alpha_end = self.angle_connections[1]
        
    def init_primitives( self ):
        
        self.capacitor = Connector_large_capacitor(DPoint(0,0), self. width, self.L_sep, self.gap_1,
                                                self.gap_2, self.spacing, self.d_1, self.d_2, self.N)
        self.primitives["capacitor"] = self.capacitor
             
#Purcell filter cavity        
        for i in range(self.N_wigling):
            name = "coilf" + str(i+1)
            setattr( self, name, Coil_type_2( self.Z0, self.capacitor.end + DPoint( (i)*(4*self.r), 0 ),
                                            self.L_wigling, self.r, self.L_wigling ) )
            self.primitives[name] = getattr( self, name )
            
        
        self.arcf = CPW_arc_2( self.Z0, self.primitives["coilf" + str(self.N_wigling)].end, -self.L_wigling - self.r, pi/2, 0 )
        self.primitives["arcf"] = self.arcf
        
        self.lengthf = CPW(self.Z0.width, self.Z0.gap, self.arcf.end, self.arcf.end + DPoint(self.L2_f, 0))
        self.primitives["lengthf"] = self.lengthf
        
# Now we create the connection between the Purcell-filter cavity and the resonator 
        
        self.L_connect_pad = (self.L_connector - self.Sep_connector)/2 + self.Z0.gap
        
        self.connector1 = CPW(self.Z0.width, self.Z0.gap, self.arcf.end + DPoint(self.Z0.width/2, -self.Z0.width/2),
                            self.arcf.end + DPoint(self.Z0.width/2, -self.Z0.width/2) + DPoint(0, -self.L_connect_pad))
        self.primitives["connector1"] = self.connector1
        
        
# Resonator cavity to which we couple the qubit for measurement 
        
        self.p0 = self.capacitor.end - DPoint(0, 2*self.r + self.Z0.width + 2*self.Z0.gap +2*self.L_wigling + self.L_connector)
        
        for i in range(self.N_wigling):
            name = "coilr" + str(i+1)
            setattr( self, name, Coil_type_2( self.Z0, self.p0 + DPoint( (i)*(4*self.r), 0 ), self.L_wigling, self.r,
                                             self.L_wigling, initial_direction = "up" ) )
            self.primitives[name] = getattr( self, name )
        
        self.arcr = CPW_arc_2( self.Z0, self.primitives["coilr" + str(self.N_wigling)].end, -self.L_wigling - self.r, -pi/2, 0 )
        self.primitives["arcr"] = self.arcr
        
        self.lengthr = CPW(self.Z0.width, self.Z0.gap, self.arcr.end, self.arcr.end + DPoint(self.L2_r, 0))
        self.primitives["lengthr"] = self.lengthr
        
        self.connector2 = CPW(self.Z0.width, self.Z0.gap, self.arcr.end + DPoint(self.Z0.width/2, self.Z0.width/2),
                             self.arcr.end + DPoint(self.Z0.width/2, self.Z0.width/2) + DPoint(0, self.L_connect_pad))
        self.primitives["connector2"] = self.connector2
        
        self.void = CPW(0, self.Z0.gap + self.Z0.width/2, self.connector1.end, self.connector2.end)
        self.primitives["void"] = self.void
        
        self.QbCPW = CPW(self.Z0.width, self.Z0.gap, self.primitives["coilr" + str(1)].start,
                         self.primitives["coilr" + str(1)].start - DPoint(0, self.L_qubit_coupling))
        self.primitives["QbCPW"] = self.QbCPW
        
        self.cop_open_end = CPW( 0, self.Z0.b/2, self.QbCPW.end, self.QbCPW.end - DPoint(0,self.Z0.b) )
        self.primitives["cop_open_end"] = self.cop_open_end        
        
        self.connections = [DPoint(0,0), self.QbCPW.end]
        self.angle_connections = [0,self.QbCPW.alpha_end]