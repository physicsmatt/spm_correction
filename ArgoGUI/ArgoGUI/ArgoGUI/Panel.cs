using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Linq;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.Data;
using FreeImageAPI;
using System.Diagnostics;
using System.Threading;
using System.ComponentModel;
using System.IO;
using System.Windows.Media.Imaging;
namespace ArgoGUI {
    public partial class Panel : Form {
        public Panel() {
            InitializeComponent();
        }

        private void Panel_Load(object sender, System.EventArgs e) {
            // throw new System.NotImplementedException();
        }

        /**
         * Outputs .tiff file to be read into Argo Solution as input parameter
         * @param layer the channel corrected
         * @return true if successful
         **/
        private bool outputBase(int layer) {
            float[] array;
            int width = base_dimensions[0];
            int height = base_dimensions[1];
            float value;

            if (layer == -1) {
                return false;
            }

            string image_name = (string)dropdown_imageLayer.Items[layer];
            array = new float[width * height];
            Array.Copy(base_data, width * height * layer, array, 0, width * height);

            FREE_IMAGE_FORMAT format = FREE_IMAGE_FORMAT.FIF_TIFF;
            FREE_IMAGE_TYPE type = FREE_IMAGE_TYPE.FIT_FLOAT;
            FREE_IMAGE_SAVE_FLAGS save_flags = FREE_IMAGE_SAVE_FLAGS.TIFF_NONE;
            FreeImageBitmap image = new FreeImageBitmap(width, height, type);

            for (int j = 0; j < height; ++j) {
                Scanline<float> scanline = image.GetScanline<float>(j);
                for (int k = 0; k < width; ++k) {
                    value = array[width * j + k];
                    scanline.SetValue(value, k);
                }
            }
            image.RotateFlip(base_orientation);
            image.Save("./Workspace/" + image_name + "_base.tiff", format, save_flags);
            return true;
        }

        /**
         * Outputs displayable image of image channel to be corrected.
         * @param layer the channel corrected
         * @return true if successful
        **/
        private void outputDisplayableBase(int layer) {
            float[] array;
            int width = base_dimensions[0];
            int height = base_dimensions[1];
            float value, value16bpp;

            string image_name = (string)dropdown_imageLayer.Items[layer];
            array = new float[width * height];
            Array.Copy(base_data, width * height * layer, array, 0, width * height);

            FREE_IMAGE_FORMAT format = FREE_IMAGE_FORMAT.FIF_TIFF;
            FREE_IMAGE_TYPE type = FREE_IMAGE_TYPE.FIT_FLOAT;
            FREE_IMAGE_SAVE_FLAGS save_flags = FREE_IMAGE_SAVE_FLAGS.TIFF_NONE;
            FreeImageBitmap image = new FreeImageBitmap(width, height, type);

            float maxValue = Enumerable.Range(0, array.Length).Max(m => array[m]);
            float minValue = Enumerable.Range(0, array.Length).Min(m => array[m]);
            for (int j = 0; j < height; ++j) {
                Scanline<float> scanline = image.GetScanline<float>(j);
                for (int k = 0; k < width; ++k) {
                    value = array[width * j + k];
                    value16bpp = (1.0f / (maxValue - minValue) * (value - minValue));
                    scanline.SetValue(value16bpp, k);
                }
            }
            image.RotateFlip(base_orientation);
            image.Save("./Workspace/" + image_name + "_base_displayable.tiff", format, save_flags);
        }

        /**
         * Outputs .tiff file to be read into Argo Solution as input parameter
         * @param layer the channel corrected
         * @true if successful
         **/
        private bool outputSliver(int layer) {
            float[] array;
            int width = sliver_dimensions[0];
            int height = sliver_dimensions[1];
            float value;

            if (layer == -1) {
                return false;
            }

            string image_name = (string)dropdown_imageLayer.Items[layer];
            array = new float[width * height];
            Array.Copy(sliver_data, width * height * layer, array, 0, width * height);

            FREE_IMAGE_FORMAT format = FREE_IMAGE_FORMAT.FIF_TIFF;
            FREE_IMAGE_TYPE type = FREE_IMAGE_TYPE.FIT_FLOAT;
            FREE_IMAGE_SAVE_FLAGS save_flags = FREE_IMAGE_SAVE_FLAGS.TIFF_NONE;
            FreeImageBitmap image = new FreeImageBitmap(width, height, type);

            for (int j = 0; j < height; ++j) {
                Scanline<float> scanline = image.GetScanline<float>(j);
                for (int k = 0; k < width; ++k) {
                    value = array[width * j + k];
                    scanline.SetValue(value, k);
                }
            }
            image.RotateFlip(sliver_orientation);
            image.Save("./Workspace/" + image_name + "_sliver.tiff", format, save_flags);
            return true;
        }

        /**
         * Outputs displayable image of image channel to be corrected.
         * @param layer the channel corrected
         * @return true if successful
        **/
        private void outputDisplayableSliver(int layer) {
            float[] array;
            int width = sliver_dimensions[0];
            int height = sliver_dimensions[1];
            float value, value16bpp;

            string image_name = (string)dropdown_imageLayer.Items[layer];
            array = new float[width * height];
            Array.Copy(sliver_data, width * height * layer, array, 0, width * height);

            FREE_IMAGE_FORMAT format = FREE_IMAGE_FORMAT.FIF_TIFF;
            FREE_IMAGE_TYPE type = FREE_IMAGE_TYPE.FIT_FLOAT;
            FREE_IMAGE_SAVE_FLAGS save_flags = FREE_IMAGE_SAVE_FLAGS.TIFF_NONE;
            FreeImageBitmap image = new FreeImageBitmap(width, height, type);

            float maxValue = Enumerable.Range(0, array.Length).Max(m => array[m]);
            float minValue = Enumerable.Range(0, array.Length).Min(m => array[m]);
            for (int j = 0; j < height; ++j) {
                Scanline<float> scanline = image.GetScanline<float>(j);
                for (int k = 0; k < width; ++k) {
                    value = array[width * j + k];
                    value16bpp = (1.0f / (maxValue - minValue) * (value - minValue));
                    scanline.SetValue(value16bpp, k);
                }
            }
            image.RotateFlip(sliver_orientation);
            image.Save("./Workspace/" + image_name + "_sliver_displayable.tiff", format, save_flags);
        }

       /**
        * Output .ibw file of the base image with selected layers corrected
       **/
        private void outputIBWBase()
        {
            int dimMult = base_dimensions[0] * base_dimensions[1] * base_dimensions[2];
            NativeCalls.SetRawData(base_name + ".ibw", (uint)dimMult, base_data);
        }




        //GH 2016, need to update base_data when the Argo runs
        /**
         * Updates base_data with the corrected float arrays from the Argo Solution for each layer
         * @param image_name name of the .tiff file in Workspace corrected by Argo Solution
         * @param layer channel of the image that was corrected and is now updated
         * @return true if successfull
        **/
        private bool updateBaseData(string image_name, int layer)
        {
            if (layer == -1)
            {
                return false;
            }
            float[] pixel_data;
            FIBITMAP image = FreeImage.Load(FREE_IMAGE_FORMAT.FIF_TIFF, image_name, FREE_IMAGE_LOAD_FLAGS.DEFAULT);
            FREE_IMAGE_TYPE type = FreeImage.GetImageType(image);
            uint width = FreeImage.GetWidth(image);
            uint height = FreeImage.GetHeight(image);
            uint pixels = width * height;
            float[] data = new float[width * height * base_dimensions[2]];
            unsafe
            {
                pixel_data = new float[pixels];
                for (int y = 0; y < height; ++y)
                {
                    float* bits = (float*)FreeImage.GetScanLine(image, y);
                    for (uint x = 0; x < width; ++x)
                    {
                        pixel_data[x + width * y] = bits[x];
                    }
                }
            }

            pixel_data.CopyTo(base_data, width * height * layer);
  
            return true;
        }

        /**
         * Moving and Deleting files after Argo Solution completion
         * @param image_name name of the channel corrected
         * @param layer layer corrected
         **/
        private void manageFiles( string image_name, int layer) {
         
            try {
                File.Move( "./base.correctedslowZ.tiff", "./Workspace/" + base_name + "/" + image_name + ".tif" ); //  + ".base.corrected.tiff" );
            }
            catch ( Exception ) {

            }
            try {
                File.Move( "./sliver.correctedslowZ.tiff", "./Workspace/" + sliver_name + "/" + image_name + ".tif" ); // +  ".sliver.corrected.tiff" );
            }
            catch ( Exception ) {

            }

            try {
                File.Delete( "./base_displayable.correctedslowZ.tiff" );
            }
            catch ( Exception ) {

            }
            try {
                File.Delete( "./sliver_displayable.correctedslowZ.tiff" );
            }
            catch ( Exception ) {

            }

            try {
                File.Move( "./base.correctedfastZ.tiff", "./Workspace/" + base_name + "/" + image_name + ".tif" ); // + ".base.corrected.tiff" );
            }
            catch ( Exception e ) {
                
            }
            try {
                File.Move( "./sliver.correctedfastZ.tiff", "./Workspace/" + sliver_name + "/" + image_name + ".tif" ); // + ".sliver.corrected.tiff" );
            }
            catch ( Exception ) {

            }

            try {
               File.Delete( "./base_displayable.correctedfastZ.tiff" );
            }
            catch ( Exception ) {

            }
            try {
                File.Delete( "./sliver_displayable.correctedfastZ.tiff" );
            }
            catch ( Exception ) {

            }

            try {
                File.Delete( "./output.txt" );

            }
            catch ( Exception ) {

            }


            try {
                File.Delete( "./Workspace/" + image_name + "_base.tiff" );
            }
            catch ( Exception ) {

            }
            try {
                File.Delete( "./Workspace/" + image_name + "_sliver.tiff" );
            }
            catch ( Exception ) {

            }
            //GH 2016    
           updateBaseData("./Workspace/" + base_name + "/" + image_name + ".tif", layer);
        }

        /**
         * Get layer to be corrected from order list
         * @param index index in order that corresponds to the layer to be corrected
         * @return the layer to be corrected
         **/
        private int getLayerForOrder( int index ) {
            return order[ index ];
        }

        /**
         * Constructs parameter string to be read into Argo Solution
         * @param image_name names of base and sliver .tiff files to be read into Argo Solution
         * @param mode designated whether working in slow-z or fast-z correction
         * @return string as argument for Argo Solution
         **/
        private string constructParameterString( string image_name, int mode ) {
            double A0Max = Convert.ToDouble( nud_A0Max.Value );
            double B0Max = Convert.ToDouble( nud_B0Max.Value );

            double a1Max = Convert.ToDouble( nud_a1Max.Value );
            double b1Max = Convert.ToDouble( nud_b1Max.Value );

            double a2Multiplier = Convert.ToDouble( nud_a2Multiplier.Value );
            double b2Multiplier = Convert.ToDouble( nud_b2Multiplier.Value );

            double a3Multiplier = Convert.ToDouble( nud_a3Multiplier.Value );
            double b3Multiplier = Convert.ToDouble( nud_b3Multilplier.Value );

            int precision = Convert.ToInt32( nud_precision.Value );
            int blockSize = Convert.ToInt32( nud_blockSize.Value );

            double growth = Convert.ToDouble( nud_growth.Value );
            double contraction = Convert.ToDouble( nud_contraction.Value );
            double reflection = Convert.ToDouble( nud_reflection.Value );
            int iterations = Convert.ToInt32( nud_iterations.Value );

            int simplex_mode;
            if(button_slowZ.Checked)    {
                simplex_mode = 0;
            }
            else   {
                simplex_mode = 1;
            }

            string commands = 
                "-i " + "./Workspace/" + image_name + "_base.tiff" + " " +
                "-s " + "./Workspace/" + image_name + "_sliver.tiff" + " " +
                "-m " + mode + " " +
                "-z " + simplex_mode + " " +
                "-0 " + A0Max + " " +
                "-1 " + B0Max + " " +
                "-2 " + a1Max + " " +
                "-3 " + b1Max + " " +
                "-4 " + a2Multiplier + " " +
                "-5 " + b2Multiplier + " " +
                "-6 " + a3Multiplier + " " +
                "-7 " + b3Multiplier + " " +
                "-p " + precision + " " +
                "-b " + blockSize + " " +
                "-g " + growth + " " +
                "-c " + contraction + " " +
                "-r " + reflection + " " +
                "-t " + iterations;

            return commands;
        }

        /**
         * Confirms that buttons have been checked, and sets the order list
         * @return false if fast-z or slow-z button has not been clicked, or if no layers are selected
         **/
        private bool validateData() {
            if ( !( button_fastZ.Checked || button_slowZ.Checked ) ) {
                return false;
            }
            if ( !( !button_baseload.Enabled && !button_sliverload.Enabled ) ) {
                return false;
            }

            bool selection = false;
            order = Enumerable.Repeat<int>( -1, checkedList_layers.Items.Count ).ToArray<int>();
            for ( int i = 0; i < checkedList_layers.Items.Count; ++i ) {
                int dash_index = ( ( String ) ( checkedList_layers.Items[ i ] ) ).LastIndexOf( " - " );
                if ( dash_index == -1 ) {
                    continue;
                }
                selection = true;
                int index = Convert.ToInt32( ( ( String ) ( checkedList_layers.Items[ i ] ) ).Substring( dash_index + 3 ) ) - 1;
                order[ index ] = i;
            }
            if ( !selection ) {
                order = null;
                return false;
            }
            return true;
        }

        //COME BACK TO THIS
        /**
         * Sets up workers to begin Argo Solution and creates directories for files
         * @param sender 
         * @param e
         **/
        private void button_correct_MouseClick( object sender, System.Windows.Forms.MouseEventArgs e ) {
            if ( e.Button == MouseButtons.Left ) {
                if ( button_correct.Text.Equals("Cancel") ) {
                    try {
                        worker.CancelAsync();
                    }
                    catch ( Exception ) {

                    }

                    button_correct.Text = "Correct!";
                    button_correct.Update();
                }
                else if ( validateData() ) {
                    System.IO.Directory.CreateDirectory( ".\\Workspace" );
                    System.IO.Directory.CreateDirectory( ".\\Workspace\\" + base_name + "\\" );
                    System.IO.Directory.CreateDirectory( ".\\Workspace\\" + sliver_name + "\\" );

                    progressbar.Value = 1;

                    groupbox_correctionMode.Enabled = false;

                    button_flipSliver.Enabled = false;
                    button_rotateSliver.Enabled = false;
                    button_flipBase.Enabled = false;
                    button_rotateBase.Enabled = false;
                    button_cleardata.Enabled = false;
                    button_defaultParameters.Enabled = false;

                    groupbox_simplexParameters.Enabled = false;
                    groupbox_gridSearchParameters.Enabled = false;
                    checkedList_layers.Enabled = false;

                    button_correct.Text = "Cancel";
                    button_correct.Update();
                    worker = new BackgroundWorker();
                    worker.WorkerReportsProgress = true;
                    worker.WorkerSupportsCancellation = true;
                    // define the event handlers
                    worker.DoWork += new DoWorkEventHandler( ExecuteSync );
                    worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler( ExecuteSyncCompleted );
                    worker.ProgressChanged += new ProgressChangedEventHandler( worker_ProgressChanged );
                    // starts the background worker
                    worker.RunWorkerAsync();
                }
                else {
                    Console.WriteLine( "Data not valid!" );
                }
            }
        }

        /**
         * Updates progress bar at bottom of form
         * @param sender
         * @param e
         **/
        private void worker_ProgressChanged( object sender, ProgressChangedEventArgs e ) {
            if ( progressbar.Value + e.ProgressPercentage >= 100 ) {
                progressbar.Value = 100;
            }
            else {
                progressbar.Value += e.ProgressPercentage;
            }
        }

        /**
         * Updates progress bar for when solution is complete, and re-enables boxes in form 
         *@param sender
         *@param args
         **/
        private void ExecuteSyncCompleted( object sender, RunWorkerCompletedEventArgs args ) {
            if ( args.Cancelled ) {
                progressbar.Value = 0;
            }
            else {
                progressbar.Value = 100;

            }
            groupbox_correctionMode.Enabled = true;

            button_flipSliver.Enabled = true;
            button_rotateSliver.Enabled = true;
            button_flipBase.Enabled = true;
            button_rotateBase.Enabled = true;
            button_cleardata.Enabled = true;
            button_defaultParameters.Enabled = true;

            groupbox_simplexParameters.Enabled = true;
            groupbox_gridSearchParameters.Enabled = true;
            checkedList_layers.Enabled = true;

            button_correct.Text = "Correct!";
            button_correct.Update();
        }

        /**
         * Begins executing Argo Solution
         * @param sender
         * @param args
         **/
        private void ExecuteSync( object sender, DoWorkEventArgs args ) {
            BackgroundWorker worker = sender as BackgroundWorker;
            for ( int i = 0; i < order.Length; ++i ) {
                if ( ( worker.CancellationPending == true ) ) {
                    args.Cancel = true;
                    break;
                }
                else {

                    int layer = getLayerForOrder( i );
                    if ( layer == -1 ) {
                        continue;
                    }
                    string image_name = ( string ) dropdown_imageLayer.Items[ layer ];
                    outputBase( layer ); 
                    outputSliver( layer );
                    
                    if (i == 0)
                    {
                        string command_args = constructParameterString(image_name, 0);  // Run fully
                        ExecuteArgoSync(command_args, worker, args);
                    }
                    else
                    {
                        string command_args = constructParameterString(image_name, 2);  // Run simplex only
                        ExecuteArgoSync(command_args, worker, args);
                    }


                    //GH 2016 added layer so we can use updateBaseData
                    manageFiles( image_name,layer);
                  //  outputDisplayableBase(layer);//GH, to be displayable
                    worker.ReportProgress( 100 / layer_order );

                }
            }
            worker.ReportProgress( 100 );
            //GH 2016
            outputIBWBase();
            try{
                File.Move(base_name + ".ibw", "./Workspace/" + base_name + ".ibw");
            }
            catch (Exception){

            }
            try{
                File.Delete(base_name + ".ibw");
            }
            catch (Exception){

            }
        }

        /** Executes Argo
         * @param args
         * @param worker
         * @param workargs
         **/
        private void ExecuteArgoSync( string args, BackgroundWorker worker, DoWorkEventArgs workargs ) {
            ProcessStartInfo startInfo = new ProcessStartInfo();
            startInfo.CreateNoWindow = false;
            startInfo.UseShellExecute = false;
            startInfo.FileName = "./Argo.exe";
            startInfo.WindowStyle = ProcessWindowStyle.Hidden;
            startInfo.Arguments = ( string ) args;

            try {
                // Start the process with the info we specified.
                // Call WaitForExit and then the using statement will close.
                using ( Process exeProcess = Process.Start( startInfo ) ) {
                    while ( !exeProcess.HasExited ) {
                        System.Threading.Thread.Sleep( 2000 );
                        if ( ( worker.CancellationPending == true ) ) {
                            workargs.Cancel = true;
                            exeProcess.Kill();
                            exeProcess.Dispose();
                            break;
                        }
                    }
                    exeProcess.Kill();
                    exeProcess.Dispose();
                }
            }
            catch ( Exception err ) {
                Console.WriteLine( err.StackTrace );
            }
        }

        /**
         * Uploads base image
         * @param sender
         * @param e
         **/
        private void button_baseload_Click( object sender, System.EventArgs e ) {
            // Create an instance of the open file dialog box.
            OpenFileDialog fileDialog = new OpenFileDialog();

            // Set filter options and filter index.
            fileDialog.Filter = "IBW Files (.ibw)|*.ibw";
            fileDialog.FilterIndex = 1;

            fileDialog.Multiselect = false;

            // Call the ShowDialog method to show the dialog box.
            DialogResult userClickedOK = fileDialog.ShowDialog( this );

            // Process input if the user clicked OK.
            if ( userClickedOK.Equals( DialogResult.OK ) ) {

                //GH 2016 added so we have the full file path, and then copy the file so we have a copy in the solution directory
                File.Copy(fileDialog.FileName, fileDialog.SafeFileName, true);

                // Open the selected file to read.
                if ( readBaseIBWFile( fileDialog.FileName ) ) {
                    base_name = fileDialog.SafeFileName;
                    base_name = base_name.Replace( ".ibw", "" );
                    Button button = ( Button ) sender;
                    button.Enabled = false;
                    button_flipBase.Enabled = true;
                    button_rotateBase.Enabled = true;

                    dropdown_imageLayer.SelectedIndex = dropdown_imageLayer.SelectedIndex != -1 ? dropdown_imageLayer.SelectedIndex : 0;
                    updateImages( dropdown_imageLayer.SelectedIndex );
                }
                else {
                    MessageBox.Show( "Base file conflicts with sliver file. Please choose another base file or clear files.", "Base selection error!",
                            MessageBoxButtons.OK, MessageBoxIcon.Asterisk );
                }
            }
        }

        /**
         * Uploads base image
         * @param sender
         * @param e
         **/
        private void button_sliverload_Click( object sender, EventArgs e ) {
            // Create an instance of the open file dialog box.
            OpenFileDialog fileDialog = new OpenFileDialog();

            // Set filter options and filter index.
            fileDialog.Filter = "IBW Files (.ibw)|*.ibw";
            fileDialog.FilterIndex = 1;

            fileDialog.Multiselect = false;

            // Call the ShowDialog method to show the dialog box.
            DialogResult userClickedOK = fileDialog.ShowDialog( this );

            // Process input if the user clicked OK.
            if ( userClickedOK.Equals( DialogResult.OK ) ) {
                // Open the selected file to read.
                if ( readSliverIBWFile( fileDialog.FileName ) ) {
                    sliver_name = fileDialog.SafeFileName;
                    sliver_name = sliver_name.Replace( ".ibw", "" );
                    Button button = ( Button ) sender;
                    button.Enabled = false;
                    button_flipSliver.Enabled = true;
                    button_rotateSliver.Enabled = true;

                    dropdown_imageLayer.SelectedIndex = dropdown_imageLayer.SelectedIndex != -1 ? dropdown_imageLayer.SelectedIndex : 0;
                    updateImages( dropdown_imageLayer.SelectedIndex );
                }
                else {
                    MessageBox.Show( "Sliver file conflicts with base file. Please choose another sliver file or clear files.", "Sliver selection error!",
                            MessageBoxButtons.OK, MessageBoxIcon.Asterisk );
                }
            }
        }

        /**
         * Reads the ibw. base file and extracts the raw data
         * @param file_name file path of the base .ibw file
         * @return true if successful
         **/
        private bool readBaseIBWFile(string file_name) {
            IntPtr ptr;
            int num_layers = NativeCalls.GetNumberOfLayers( file_name );

            if ( !button_sliverload.Enabled ) {
                if ( sliver_dimensions[ 2 ] != num_layers ) {
                    return false;
                }
                for ( int i = 0; i < num_layers; ++i ) {
                    ptr = NativeCalls.GetLayerName( file_name, i );
                    string name = Marshal.PtrToStringAnsi( ptr );
                    NativeCalls.ReleaseMemory( ptr );
                    if ( dropdown_imageLayer.Items.IndexOf( name ) != i ) {
                        return false;
                    }
                }

                base_dimensions = new int[ 3 ];
                ptr = NativeCalls.GetDimensions( file_name );
                Marshal.Copy( ptr, base_dimensions, 0, 3 );
                NativeCalls.ReleaseMemory( ptr );

                if ( base_dimensions[ 0 ] < sliver_dimensions[ 0 ] || base_dimensions[ 1 ] < sliver_dimensions[ 1 ] ) {
                    base_dimensions = null;
                    return false;
                }
                if ( ! ( base_dimensions[ 0 ] == sliver_dimensions[ 0 ] ^ base_dimensions[ 1 ] == sliver_dimensions[ 1 ] ) ) {
                    base_dimensions = null;
                    return false;
                }

                int xyMult = base_dimensions[ 0 ] * base_dimensions[ 1 ];
                int dimMult = xyMult * base_dimensions[ 2 ];
                base_data = new float[ dimMult ];
                ptr = NativeCalls.GetRawData( file_name, ( uint ) dimMult );
                Marshal.Copy( ptr, base_data, 0, dimMult );
                NativeCalls.ReleaseMemory( ptr );
                return true;
            }
            else {
                dropdown_imageLayer.Items.Clear();
                checkedList_layers.Items.Clear();
                for ( int i = 0; i < num_layers; ++i ) {
                    ptr = NativeCalls.GetLayerName( file_name, i );
                    string name = Marshal.PtrToStringAnsi( ptr );
                    NativeCalls.ReleaseMemory( ptr );
                    dropdown_imageLayer.Items.Insert( i, name );
                    checkedList_layers.Items.Insert( i, name );
                }

                base_dimensions = new int[ 3 ];
                ptr = NativeCalls.GetDimensions( file_name );
                Marshal.Copy( ptr, base_dimensions, 0, 3 );
                NativeCalls.ReleaseMemory( ptr );
                int xyMult = base_dimensions[ 0 ] * base_dimensions[ 1 ];
                int dimMult = xyMult * base_dimensions[ 2 ];
                base_data = new float[ dimMult ];
                ptr = NativeCalls.GetRawData( file_name, ( uint ) dimMult );
                Marshal.Copy( ptr, base_data, 0, dimMult );
                NativeCalls.ReleaseMemory( ptr );
                return true;
            }
        }

        /**
         * Reads the ibw. base file and extracts the raw data
         * @param file_name file path of the base .ibw file
         * @return true if successful
         **/
        private bool readSliverIBWFile( string file_name ) {
            IntPtr ptr;
            int num_layers = NativeCalls.GetNumberOfLayers( file_name );

            if ( !button_baseload.Enabled ) {
                if ( base_dimensions[ 2 ] != num_layers ) {
                    return false;
                }
                for ( int i = 0; i < num_layers; ++i ) {
                    ptr = NativeCalls.GetLayerName( file_name, i );
                    string name = Marshal.PtrToStringAnsi( ptr );
                    NativeCalls.ReleaseMemory( ptr );
                    if ( dropdown_imageLayer.Items.IndexOf( name ) != i ) {
                        return false;
                    }
                }

                sliver_dimensions = new int[ 3 ];
                ptr = NativeCalls.GetDimensions( file_name );
                Marshal.Copy( ptr, sliver_dimensions, 0, 3 );
                NativeCalls.ReleaseMemory( ptr );

                if ( base_dimensions[ 0 ] < sliver_dimensions[ 0 ] || base_dimensions[ 1 ] < sliver_dimensions[ 1 ] ) {
                    sliver_dimensions = null;
                    return false;
                }
                if ( ! ( base_dimensions[ 0 ] == sliver_dimensions[ 0 ] ^ base_dimensions[ 1 ] == sliver_dimensions[ 1 ] ) ) {
                    sliver_dimensions = null;
                    return false;
                }

                int xyMult = sliver_dimensions[ 0 ] * sliver_dimensions[ 1 ];
                int dimMult = xyMult * sliver_dimensions[ 2 ];
                sliver_data = new float[ dimMult ];
                ptr = NativeCalls.GetRawData( file_name, ( uint ) dimMult );
                Marshal.Copy( ptr, sliver_data, 0, dimMult );
                NativeCalls.ReleaseMemory( ptr );
                return true;
            }
            else {
                dropdown_imageLayer.Items.Clear();
                checkedList_layers.Items.Clear();
                for ( int i = 0; i < num_layers; ++i ) {
                    ptr = NativeCalls.GetLayerName( file_name, i );
                    string name = Marshal.PtrToStringAnsi( ptr );
                    NativeCalls.ReleaseMemory( ptr );
                    dropdown_imageLayer.Items.Insert( i, name );
                    checkedList_layers.Items.Insert( i, name );
                }

                sliver_dimensions = new int[ 3 ];
                ptr = NativeCalls.GetDimensions( file_name );
                Marshal.Copy( ptr, sliver_dimensions, 0, 3 );
                NativeCalls.ReleaseMemory( ptr );
                int xyMult = sliver_dimensions[ 0 ] * sliver_dimensions[ 1 ];
                int dimMult = xyMult * sliver_dimensions[ 2 ];
                sliver_data = new float[ dimMult ];
                ptr = NativeCalls.GetRawData( file_name, ( uint ) dimMult );
                Marshal.Copy( ptr, sliver_data, 0, dimMult );
                NativeCalls.ReleaseMemory( ptr );
                return true;
            }
        }

        /**
         */
        private void updateImages( int index ) {
            float[] array;
            float maxValue, minValue, value;
            int value16bpp;
            Bitmap bitmap;
            Color color;
            Image old_bitmap;

            if ( base_dimensions != null ) {
                int width = base_dimensions[ 0 ];
                int height = base_dimensions[ 1 ];

                array = new float[ width * height ];
                Array.Copy( base_data, width * height * index, array, 0, width * height );

                FreeImageBitmap image_bitmap = new FreeImageBitmap( new Bitmap(width, height), width, height );
                maxValue = Enumerable.Range( 0, array.Length ).Max( m => array[ m ] );
                minValue = Enumerable.Range( 0, array.Length ).Min( m => array[ m ] );
                for ( int j = 0; j < height; ++j ) {
                    Scanline<int> scanline = image_bitmap.GetScanline<int>( j );

                    for ( int k = 0; k < width; ++k ) {
                        value = array[ width * j + k ];
                        value16bpp = (int) ( 255.0f / ( maxValue - minValue ) * ( value - minValue ) );
                        color = Color.FromArgb( value16bpp, value16bpp, value16bpp );
                        scanline.SetValue( color.ToArgb(), k );
                    }
                }

                image_bitmap.RotateFlip( base_orientation );
                bitmap = image_bitmap.ToBitmap();
                image_bitmap.Dispose();

                old_bitmap = image_base.Image;
                image_base.Image = bitmap;
                if ( old_bitmap != null )
                    old_bitmap.Dispose();
            }



            if ( sliver_dimensions != null ) {
                int width = sliver_dimensions[ 0 ];
                int height = sliver_dimensions[ 1 ];

                array = new float[ width * height ];
                Array.Copy( sliver_data, width * height * index, array, 0, width * height );
                FreeImageBitmap image_bitmap = new FreeImageBitmap( new Bitmap(width, height), width, height );
                maxValue = Enumerable.Range( 0, array.Length ).Max( m => array[ m ] );
                minValue = Enumerable.Range( 0, array.Length ).Min( m => array[ m ] );
                for ( int j = 0; j < height; ++j ) {
                    Scanline<int> scanline = image_bitmap.GetScanline<int>( j );

                    for ( int k = 0; k < width; ++k ) {
                        value = array[ width * j + k ];
                        value16bpp = (int) ( 255.0f / ( maxValue - minValue ) * ( value - minValue ) );
                        color = Color.FromArgb( value16bpp, value16bpp, value16bpp );
                        scanline.SetValue( color.ToArgb(), k );
                    }
                }

                image_bitmap.RotateFlip( sliver_orientation );
                bitmap = image_bitmap.ToBitmap();
                image_bitmap.Dispose();

                old_bitmap = image_sliver.Image;
                image_sliver.Image = bitmap;
                if ( old_bitmap != null )
                    old_bitmap.Dispose();
            }
        }

        /**
         * 
         **/
        private void checkedList_layers_ItemCheck( object sender, ItemCheckEventArgs e ) {
            int index = e.Index;
            if ( e.NewValue == CheckState.Checked ) {
                // Just checked item
                layer_order++;
                checkedList_layers.Items[ index ] = checkedList_layers.Items[ index ] + " - " + layer_order;
            }
            else {
                int dash_index;
                int img1_index;
                int img2_index;

                dash_index = ( ( String ) ( checkedList_layers.Items[ index ] ) ).LastIndexOf( " - " );
                img1_index = Convert.ToInt32( ( ( String ) ( checkedList_layers.Items[ index ] ) ).Substring( dash_index + 3 ) );
                for ( int i = 0; i < checkedList_layers.Items.Count; ++i ) {
                    dash_index = ( ( String ) ( checkedList_layers.Items[ i ] ) ).LastIndexOf( " - " );
                    if ( dash_index == -1 ) {
                        continue;
                    }
                    img2_index = Convert.ToInt32( ( ( String ) ( checkedList_layers.Items[ i ] ) ).Substring( dash_index + 3 ) );
                    if ( img2_index < img1_index ) {
                        continue;
                    }
                    checkedList_layers.Items[ i ] = ( ( String ) ( checkedList_layers.Items[ i ] ) ).Remove( dash_index ) + " - " + ( img2_index - 1 );
                }
                layer_order--;
                dash_index = ( ( String ) ( checkedList_layers.Items[ index ] ) ).LastIndexOf( " - " );
                checkedList_layers.Items[ index ] = ( ( String ) ( checkedList_layers.Items[ index ] ) ).Remove( dash_index );
            }
        }

        /**
         * Check orientation of the base image shown in the window, then updates if button is clicked
         * @param sender
         * @param e
         **/
        private void button_flipBase_Click( object sender, EventArgs e ) {
            switch ( base_orientation ) {
                case RotateFlipType.RotateNoneFlipNone :
                    base_orientation = RotateFlipType.RotateNoneFlipX;
                    break;
                case RotateFlipType.Rotate90FlipNone :
                    base_orientation = RotateFlipType.Rotate90FlipX;
                    break;
                case RotateFlipType.Rotate180FlipNone :
                    base_orientation = RotateFlipType.Rotate180FlipX;
                    break;
                case RotateFlipType.Rotate270FlipNone :
                    base_orientation = RotateFlipType.Rotate270FlipX;
                    break;
                case RotateFlipType.RotateNoneFlipX :
                    base_orientation = RotateFlipType.RotateNoneFlipNone;
                    break;
                case RotateFlipType.Rotate90FlipX :
                    base_orientation = RotateFlipType.Rotate90FlipNone;
                    break;
                case RotateFlipType.Rotate180FlipX :
                    base_orientation = RotateFlipType.Rotate180FlipNone;
                    break;
                case RotateFlipType.Rotate270FlipX :
                    base_orientation = RotateFlipType.Rotate270FlipNone;
                    break;
            }

            updateImages( dropdown_imageLayer.SelectedIndex );
        }

        /**
         * Check orientation of the base image shown in the window, then updates if button is clicked
         * @param sender
         * @param e
         **/
        private void button_rotateBase_Click( object sender, EventArgs e ) {
            switch ( base_orientation ) {
                case RotateFlipType.RotateNoneFlipNone:
                    base_orientation = RotateFlipType.Rotate270FlipNone;
                    break;
                case RotateFlipType.Rotate90FlipNone:
                    base_orientation = RotateFlipType.RotateNoneFlipNone;
                    break;
                case RotateFlipType.Rotate180FlipNone:
                    base_orientation = RotateFlipType.Rotate90FlipNone;
                    break;
                case RotateFlipType.Rotate270FlipNone:
                    base_orientation = RotateFlipType.Rotate180FlipNone;
                    break;
                case RotateFlipType.RotateNoneFlipX:
                    base_orientation = RotateFlipType.Rotate90FlipX;
                    break;
                case RotateFlipType.Rotate90FlipX:
                    base_orientation = RotateFlipType.Rotate180FlipX;
                    break;
                case RotateFlipType.Rotate180FlipX:
                    base_orientation = RotateFlipType.Rotate270FlipX;
                    break;
                case RotateFlipType.Rotate270FlipX:
                    base_orientation = RotateFlipType.RotateNoneFlipX;
                    break;
            }

            updateImages( dropdown_imageLayer.SelectedIndex );
        }

        /**
         * Check orientation of the sliver image shown in the window, then updates if button is clicked
         * @param sender
         * @param e
         **/
        private void button_flipSliver_Click( object sender, EventArgs e ) {
            switch ( sliver_orientation ) {
                case RotateFlipType.RotateNoneFlipNone:
                    sliver_orientation = RotateFlipType.RotateNoneFlipX;
                    break;
                case RotateFlipType.Rotate90FlipNone:
                    sliver_orientation = RotateFlipType.Rotate90FlipX;
                    break;
                case RotateFlipType.Rotate180FlipNone:
                    sliver_orientation = RotateFlipType.Rotate180FlipX;
                    break;
                case RotateFlipType.Rotate270FlipNone:
                    sliver_orientation = RotateFlipType.Rotate270FlipX;
                    break;
                case RotateFlipType.RotateNoneFlipX:
                    sliver_orientation = RotateFlipType.RotateNoneFlipNone;
                    break;
                case RotateFlipType.Rotate90FlipX:
                    sliver_orientation = RotateFlipType.Rotate90FlipNone;
                    break;
                case RotateFlipType.Rotate180FlipX:
                    sliver_orientation = RotateFlipType.Rotate180FlipNone;
                    break;
                case RotateFlipType.Rotate270FlipX:
                    sliver_orientation = RotateFlipType.Rotate270FlipNone;
                    break;
            }

            updateImages( dropdown_imageLayer.SelectedIndex );
        }

        /**
         * Check orientation of the sliver image shown in the window, then updates if button is clicked
         * @param sender
         * @param e
         **/
        private void button_rotateSliver_Click( object sender, EventArgs e ) {
            switch ( sliver_orientation ) {
                case RotateFlipType.RotateNoneFlipNone:
                    sliver_orientation = RotateFlipType.Rotate270FlipNone;
                    break;
                case RotateFlipType.Rotate90FlipNone:
                    sliver_orientation = RotateFlipType.RotateNoneFlipNone;
                    break;
                case RotateFlipType.Rotate180FlipNone:
                    sliver_orientation = RotateFlipType.Rotate90FlipNone;
                    break;
                case RotateFlipType.Rotate270FlipNone:
                    sliver_orientation = RotateFlipType.Rotate180FlipNone;
                    break;
                case RotateFlipType.RotateNoneFlipX:
                    sliver_orientation = RotateFlipType.Rotate90FlipX;
                    break;
                case RotateFlipType.Rotate90FlipX:
                    sliver_orientation = RotateFlipType.Rotate180FlipX;
                    break;
                case RotateFlipType.Rotate180FlipX:
                    sliver_orientation = RotateFlipType.Rotate270FlipX;
                    break;
                case RotateFlipType.Rotate270FlipX:
                    sliver_orientation = RotateFlipType.RotateNoneFlipX;
                    break;
            }

            updateImages( dropdown_imageLayer.SelectedIndex );
        }

        /** Updates if layer correction order is changed
         * @param sender
         * @param e
         **/
        private void dropdown_imageLayer_SelectedIndexChanged( object sender, EventArgs e ) {
            ComboBox combos = ( ComboBox ) sender;
            int index = combos.SelectedIndex;
            updateImages( index );
        }

        /** Set parameters back to default values
         * @param sender
         * @param e
         **/
        private void button_defaultParameters_Click( object sender, EventArgs e ) {
            progressbar.Value = 0;

            nud_A0Max.Value = new Decimal( 3 );
            nud_B0Max.Value = new Decimal( 3 );

            nud_a1Max.Value = new Decimal( 0.1 );
            nud_b1Max.Value = new Decimal( 0.1 );

            nud_a2Multiplier.Value = new Decimal( 0.04 );
            nud_b2Multiplier.Value = new Decimal( 0.04 );

            nud_a3Multiplier.Value = new Decimal( 0.04 );
            nud_b3Multilplier.Value = new Decimal( 0.04 );

            nud_blockSize.Value = new Decimal( 0 );
            nud_precision.Value = new Decimal( 1 );

            nud_growth.Value = new Decimal( 1.8 );
            nud_contraction.Value = new Decimal( 0.7 );
            nud_reflection.Value = new Decimal( 1.3 );
            nud_iterations.Value = new Decimal( 10000 );
        }

        /**
         * Clears data
         * @param sender
         * @param e
         **/
        private void button_cleardata_Click( object sender, EventArgs e ) {
            this.layer_order = 0;
            this.order = null;
            this.label_correctionOrder = null;

            this.progressbar.Value = 0;
            this.base_orientation = RotateFlipType.RotateNoneFlipNone;
            this.sliver_orientation = RotateFlipType.RotateNoneFlipNone;

            this.button_flipSliver.Enabled = false;
            this.button_rotateSliver.Enabled = false;
            this.button_flipBase.Enabled = false;
            this.button_rotateBase.Enabled = false;

            this.button_sliverload.Enabled = true;
            this.button_baseload.Enabled = true;

            this.dropdown_imageLayer.Items.Clear();
            this.checkedList_layers.Items.Clear();

            this.base_name = null;
            this.sliver_name = null;
            this.base_data = null;
            this.sliver_data = null;
            this.base_dimensions = null;
            this.sliver_dimensions = null;

            this.image_base.Image = null;
            this.image_sliver.Image = null;

            button_defaultParameters_Click( sender, e );
        }
    }
}