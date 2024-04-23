# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 12:18:12 2018

@author: idariash
"""

import struct

class dataN1( object ):
    
    """ Class to hold a dataN1 object that hold the data of one ray """
        
    def __new__( cls, filePath, dataInitialPosition = 0 ):
        with open( filePath, "rb" ) as fileId:
            
            fileId.seek( 0, 2 )
            
            lastFilePosition = fileId.tell()
            
            if dataInitialPosition >= lastFilePosition:
                return None
            else:
                return object.__new__( cls )
    
    
    def __init__( self, filePath, dataInitialPosition = 0 ):        
        
        with open( filePath, "rb" ) as fileId:
            
            fileId.seek( 0, 2 )
            
            lastFilePosition = fileId.tell()
            
            if dataInitialPosition >= lastFilePosition:
                return None
             
            fileId.seek( dataInitialPosition )
            self.initialPosition = fileId.tell()
            
            self.version = struct.unpack( ' h', fileId.read( 2 ) )[0]
            self.drxVersion = struct.unpack( ' h', fileId.read( 2 ) )[0]
            _ = fileId.read( 4 )
            self.initCW = struct.unpack( ' d', fileId.read( 8 ) )[0]
            self.azimuth = struct.unpack( ' f', fileId.read( 4 ) )[0]
            self.elevation = struct.unpack( ' f', fileId.read( 4 ) )[0]
            self.idVolume = struct.unpack( ' H', fileId.read( 2 ) )[0]
            self.idScan = struct.unpack( ' H', fileId.read( 2 ) )[0]
            self.idSet = struct.unpack( ' H', fileId.read( 2 ) )[0]
            self.idGroup = struct.unpack( ' H', fileId.read( 2 ) )[0]
            self.idPulse = struct.unpack( ' H', fileId.read( 2 ) )[0]
            self.scanInit = struct.unpack( ' ?', fileId.read( 1 ) )[0]
            self.scanEnd = struct.unpack( ' ?', fileId.read( 1 ) )[0]
            self.groupEnd = struct.unpack( ' ?', fileId.read( 1 ) )[0]
            self.inhibit = struct.unpack( ' ?', fileId.read( 1 ) )[0]
            self.validSamples = struct.unpack( ' H', fileId.read( 2 ) )[0]
            self.aquisitionNumber = struct.unpack( ' H', fileId.read( 2 ) )[0]
            _ = fileId.read( 2 )
            self.sequenceNumber = struct.unpack( ' I', fileId.read( 4 ) )[0]
            
            self.readTimeHigh = struct.unpack( ' Q', fileId.read( 8 ) )[0]  # seconds since epoch in ArgLocalTime or something similar
            
            self.readTimeLow = struct.unpack( ' Q', fileId.read( 8 ) )[0]  # nanoseconds since epoch in ArgLocalTime
                        
            _ = fileId.read( 64 )
            
            self.timeStampSecs = int( self.readTimeHigh ) + 3600 * 4  # In UTC (correction)
            self.timeStampNanoSecs = int( self.readTimeLow )             
            self.timeStamp = str( self.readTimeHigh ) + "%06d" % self.readTimeLow

            self.V_I = list()
            self.V_Q = list()
            self.H_I = list()
            self.H_Q = list()
    
            for _ in range( 0, self.validSamples ):
                self.V_I.append( struct.unpack( ' f' , fileId.read( 4 ) ) )
                self.V_Q.append( struct.unpack( ' f' , fileId.read( 4 ) ) )
                self.H_I.append( struct.unpack( ' f' , fileId.read( 4 ) ) )
                self.H_Q.append( struct.unpack( ' f' , fileId.read( 4 ) ) )
                
            self.finalPosition = fileId.tell() - 1 
            
    
    def dump( self, data = True ):
       
        print(" ")
        print('version =  ' + str( self.version ))
        print('drxVersion =  ' + str( self.drxVersion ))
        
        print('initCW =  ' + str( self.initCW ))
        print('azimuth =  ' + str( self.azimuth ))
        print('elevation =  ' + str( self.elevation ))
        print('idVolume =  ' + str( self.idVolume ))
        print('idScan =  ' + str( self.idScan ))
        print('idSet =  ' + str( self.idSet ))
        print('idGroup =  ' + str( self.idGroup ))
        print('idPulse =  ' + str( self.idPulse ))
        print('scanInit =  ' + str( self.scanInit ))
        print('scanEnd =  ' + str( self.scanEnd ))
        print('groupEnd =  ' + str( self.groupEnd ))
        print('inhibit =  ' + str( self.inhibit ))
        print('validSamples =  ' + str( self.validSamples ))
        print('aquisitionNumber =  ' + str( self.aquisitionNumber ))
        
        print('sequenceNumber =  ' + str( self.sequenceNumber ))
        print('readTimeHigh =  ' + str( self.readTimeHigh ))
        print('readTimeLow =  ' + str( self.readTimeLow ))
        
        print('timestampSecs =  ' + str( self.timeStampSecs ))
        print('timestampNanoSecs =  ' + str( self.timeStampNanoSecs ))
        
        print('timestamp =  ' + str( self.timeStamp ))
                
        if data:
    
            for sample in range( 0, self.validSamples ):
                print('V_I[' + str( sample ) + '] =  ' + str( self.V_I[sample] ))
                print('V_Q[' + str( sample ) + '] =  ' + str( self.V_Q[sample] ))
                print('H_I[' + str( sample ) + '] =  ' + str( self.H_I[sample] ))
                print('H_Q[' + str( sample ) + '] =  ' + str( self.H_Q[sample] ))
                
        print("dataInitialPosition =" + str( self.initialPosition )) 
        print("dataFinalPosition =" + str( self.finalPosition ))


class invapIQ( object ):
    '''
    Class to read and IQ file
    '''
   
    def __init__( self, filePath ):
        '''
        Constructor
        '''        
        
        self.filePath = filePath
        self.IQdata = list()
        
        print("Reading IQ data")
        lastIndex = 0
        while ( True ):  
            
            _data1 = dataN1( filePath, lastIndex )
            #print(str(lastIndex))
            
            if _data1 is not None:
                self.IQdata.append( _data1 )
                lastIndex = _data1.finalPosition + 1 
            else:
                break    
            
        print("Sorting by pulse ID")
        self.IQdata.sort( key = lambda x: x.timeStamp )
    
    def getScans( self ):
        """ Get the elevations available in the file """
        
        scans = list()
        
        for _dataN1 in self.IQdata:
            if _dataN1.elevation not in scans:
                scans.append( _dataN1.elevation )
        
        print(scans)
        
    def dump( self ):
        """ Dump into screen all the data read """
        for _dataN1 in self.IQdata:
            _dataN1.dump()