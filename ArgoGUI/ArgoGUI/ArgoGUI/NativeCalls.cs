
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace ArgoGUI {
    /**
     * Class used to import TestWave.dll file to be used in reading from and writing to .ibw files
     **/
    public class NativeCalls {


        
        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int DoAsylumReadTest( [MarshalAs( UnmanagedType.LPStr )] string c );

        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetNumberOfLayers( [MarshalAs( UnmanagedType.LPStr )] string c );

        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr GetLayerName( [MarshalAs( UnmanagedType.LPStr )] string c, int index );

        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr GetDimensions( [MarshalAs( UnmanagedType.LPStr )] string c );

        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern IntPtr GetRawData( [MarshalAs( UnmanagedType.LPStr )] string c, uint size );

        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ReleaseMemory( IntPtr ptr );

        [DllImport("TestWave.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int SetRawData([MarshalAs(UnmanagedType.LPStr)] string c, uint size, float[] rawData);
        
    }
}
