using System.ComponentModel;
using System.Diagnostics;
using System.Drawing;
using System.Windows.Forms;

namespace ArgoGUI
{
    partial class Panel
    {
        BackgroundWorker worker;

        private int layer_order = 0;
        private int[] order;

        private string base_name;
        private string sliver_name;
        private int[] base_dimensions;
        private int[] sliver_dimensions;
        private float[] base_data;
        private float[] sliver_data;
        private RotateFlipType base_orientation = RotateFlipType.RotateNoneFlipNone;
        private RotateFlipType sliver_orientation = RotateFlipType.RotateNoneFlipNone;

        private RadioButton button_fastZ;
        private RadioButton button_slowZ;

        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        private ComboBox dropdown_imageLayer;
        private CheckedListBox checkedList_layers;
        private GroupBox groupbox_correctionMode;
        private PictureBox image_sliver;
        private PictureBox image_base;
        private ProgressBar progressbar;

        private GroupBox groupbox_gridSearchParameters;
        private Label label_a2Multiplier;
        private Label label_a1Max;
        private Label label_A0Max;
        private Label label_a3Multiplier;
        private Label label_b3Multiplier;
        private Label label_b2Multiplier;
        private Label label_b1Max;
        private Label label_B0Max;
        private Label label_blockSize;
        private Label label_precision;
        private GroupBox groupbox_simplexParameters;
        private Label label_iterations;
        private Label label_reflection;
        private Label label_contraction;
        private Label label_growth;
        private Label label_correctionOrder;


        private Button button_correct;
        private Button button_baseload;
        private Button button_sliverload;
        private Button button_cleardata;
        private Button button_flipBase;
        private Button button_rotateBase;
        private Button button_rotateSliver;
        private Button button_flipSliver;
        private Button button_defaultParameters;


        private NumericUpDown nud_blockSize;
        private NumericUpDown nud_precision;
        private NumericUpDown nud_b3Multilplier;
        private NumericUpDown nud_b2Multiplier;
        private NumericUpDown nud_b1Max;
        private NumericUpDown nud_B0Max;
        private NumericUpDown nud_a3Multiplier;
        private NumericUpDown nud_a2Multiplier;
        private NumericUpDown nud_a1Max;
        private NumericUpDown nud_A0Max;
        private NumericUpDown nud_iterations;
        private NumericUpDown nud_reflection;
        private NumericUpDown nud_contraction;
        private NumericUpDown nud_growth;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.Windows.Forms.Label label_imageLayer;
            System.Windows.Forms.Label label_baseImage;
            System.Windows.Forms.Label label_sliverImage;
            this.dropdown_imageLayer = new System.Windows.Forms.ComboBox();
            this.groupbox_correctionMode = new System.Windows.Forms.GroupBox();
            this.button_slowZ = new System.Windows.Forms.RadioButton();
            this.button_fastZ = new System.Windows.Forms.RadioButton();
            this.progressbar = new System.Windows.Forms.ProgressBar();
            this.button_baseload = new System.Windows.Forms.Button();
            this.button_correct = new System.Windows.Forms.Button();
            this.groupbox_gridSearchParameters = new System.Windows.Forms.GroupBox();
            this.nud_blockSize = new System.Windows.Forms.NumericUpDown();
            this.nud_precision = new System.Windows.Forms.NumericUpDown();
            this.nud_b3Multilplier = new System.Windows.Forms.NumericUpDown();
            this.nud_b2Multiplier = new System.Windows.Forms.NumericUpDown();
            this.nud_b1Max = new System.Windows.Forms.NumericUpDown();
            this.nud_B0Max = new System.Windows.Forms.NumericUpDown();
            this.nud_a3Multiplier = new System.Windows.Forms.NumericUpDown();
            this.nud_a2Multiplier = new System.Windows.Forms.NumericUpDown();
            this.nud_a1Max = new System.Windows.Forms.NumericUpDown();
            this.nud_A0Max = new System.Windows.Forms.NumericUpDown();
            this.label_blockSize = new System.Windows.Forms.Label();
            this.label_precision = new System.Windows.Forms.Label();
            this.label_b3Multiplier = new System.Windows.Forms.Label();
            this.label_b2Multiplier = new System.Windows.Forms.Label();
            this.label_b1Max = new System.Windows.Forms.Label();
            this.label_B0Max = new System.Windows.Forms.Label();
            this.label_a3Multiplier = new System.Windows.Forms.Label();
            this.label_a2Multiplier = new System.Windows.Forms.Label();
            this.label_a1Max = new System.Windows.Forms.Label();
            this.label_A0Max = new System.Windows.Forms.Label();
            this.groupbox_simplexParameters = new System.Windows.Forms.GroupBox();
            this.nud_iterations = new System.Windows.Forms.NumericUpDown();
            this.label_iterations = new System.Windows.Forms.Label();
            this.nud_reflection = new System.Windows.Forms.NumericUpDown();
            this.label_reflection = new System.Windows.Forms.Label();
            this.nud_contraction = new System.Windows.Forms.NumericUpDown();
            this.label_contraction = new System.Windows.Forms.Label();
            this.nud_growth = new System.Windows.Forms.NumericUpDown();
            this.label_growth = new System.Windows.Forms.Label();
            this.button_sliverload = new System.Windows.Forms.Button();
            this.button_cleardata = new System.Windows.Forms.Button();
            this.label_correctionOrder = new System.Windows.Forms.Label();
            this.checkedList_layers = new System.Windows.Forms.CheckedListBox();
            this.button_rotateSliver = new System.Windows.Forms.Button();
            this.button_flipSliver = new System.Windows.Forms.Button();
            this.button_rotateBase = new System.Windows.Forms.Button();
            this.button_flipBase = new System.Windows.Forms.Button();
            this.image_sliver = new System.Windows.Forms.PictureBox();
            this.image_base = new System.Windows.Forms.PictureBox();
            this.button_defaultParameters = new System.Windows.Forms.Button();
            label_imageLayer = new System.Windows.Forms.Label();
            label_baseImage = new System.Windows.Forms.Label();
            label_sliverImage = new System.Windows.Forms.Label();
            this.groupbox_correctionMode.SuspendLayout();
            this.groupbox_gridSearchParameters.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_blockSize)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_precision)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_b3Multilplier)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_b2Multiplier)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_b1Max)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_B0Max)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_a3Multiplier)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_a2Multiplier)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_a1Max)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_A0Max)).BeginInit();
            this.groupbox_simplexParameters.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_iterations)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_reflection)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_contraction)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_growth)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.image_sliver)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.image_base)).BeginInit();
            this.SuspendLayout();
            // 
            // label_imageLayer
            // 
            label_imageLayer.Anchor = System.Windows.Forms.AnchorStyles.Left;
            label_imageLayer.AutoSize = true;
            label_imageLayer.Location = new System.Drawing.Point(431, 95);
            label_imageLayer.Name = "label_imageLayer";
            label_imageLayer.Size = new System.Drawing.Size(68, 13);
            label_imageLayer.TabIndex = 0;
            label_imageLayer.Text = "Image Layer:";
            // 
            // label_baseImage
            // 
            label_baseImage.AutoSize = true;
            label_baseImage.Location = new System.Drawing.Point(12, 9);
            label_baseImage.Name = "label_baseImage";
            label_baseImage.Size = new System.Drawing.Size(69, 13);
            label_baseImage.TabIndex = 0;
            label_baseImage.Text = "Base Image :";
            // 
            // label_sliverImage
            // 
            label_sliverImage.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            label_sliverImage.AutoSize = true;
            label_sliverImage.Location = new System.Drawing.Point(748, 9);
            label_sliverImage.Name = "label_sliverImage";
            label_sliverImage.Size = new System.Drawing.Size(71, 13);
            label_sliverImage.TabIndex = 0;
            label_sliverImage.Text = "Sliver Image :";
            // 
            // dropdown_imageLayer
            // 
            this.dropdown_imageLayer.Anchor = System.Windows.Forms.AnchorStyles.Left;
            this.dropdown_imageLayer.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.dropdown_imageLayer.FormattingEnabled = true;
            this.dropdown_imageLayer.ImeMode = System.Windows.Forms.ImeMode.Disable;
            this.dropdown_imageLayer.Location = new System.Drawing.Point(434, 111);
            this.dropdown_imageLayer.MaxDropDownItems = 50;
            this.dropdown_imageLayer.Name = "dropdown_imageLayer";
            this.dropdown_imageLayer.Size = new System.Drawing.Size(233, 21);
            this.dropdown_imageLayer.TabIndex = 1;
            this.dropdown_imageLayer.SelectedIndexChanged += new System.EventHandler(this.dropdown_imageLayer_SelectedIndexChanged);
            // 
            // groupbox_correctionMode
            // 
            this.groupbox_correctionMode.Anchor = System.Windows.Forms.AnchorStyles.Left;
            this.groupbox_correctionMode.Controls.Add(this.button_slowZ);
            this.groupbox_correctionMode.Controls.Add(this.button_fastZ);
            this.groupbox_correctionMode.Location = new System.Drawing.Point(434, 25);
            this.groupbox_correctionMode.Name = "groupbox_correctionMode";
            this.groupbox_correctionMode.Size = new System.Drawing.Size(200, 67);
            this.groupbox_correctionMode.TabIndex = 4;
            this.groupbox_correctionMode.TabStop = false;
            this.groupbox_correctionMode.Text = "Correction Mode";
            // 
            // button_slowZ
            // 
            this.button_slowZ.AutoSize = true;
            this.button_slowZ.Location = new System.Drawing.Point(6, 19);
            this.button_slowZ.Name = "button_slowZ";
            this.button_slowZ.Size = new System.Drawing.Size(131, 17);
            this.button_slowZ.TabIndex = 1;
            this.button_slowZ.TabStop = true;
            this.button_slowZ.Text = "Use Slow-Z Correction";
            this.button_slowZ.UseVisualStyleBackColor = true;
            // 
            // button_fastZ
            // 
            this.button_fastZ.AutoSize = true;
            this.button_fastZ.Location = new System.Drawing.Point(6, 42);
            this.button_fastZ.Name = "button_fastZ";
            this.button_fastZ.Size = new System.Drawing.Size(128, 17);
            this.button_fastZ.TabIndex = 2;
            this.button_fastZ.Text = "Use Fast-Z Correction";
            this.button_fastZ.UseVisualStyleBackColor = true;
            // 
            // progressbar
            // 
            this.progressbar.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.progressbar.Location = new System.Drawing.Point(12, 589);
            this.progressbar.Name = "progressbar";
            this.progressbar.Size = new System.Drawing.Size(717, 23);
            this.progressbar.TabIndex = 9;
            // 
            // button_baseload
            // 
            this.button_baseload.Location = new System.Drawing.Point(428, 354);
            this.button_baseload.Name = "button_baseload";
            this.button_baseload.Size = new System.Drawing.Size(90, 23);
            this.button_baseload.TabIndex = 2;
            this.button_baseload.Text = "Load Base";
            this.button_baseload.UseVisualStyleBackColor = true;
            this.button_baseload.Click += new System.EventHandler(this.button_baseload_Click);
            // 
            // button_correct
            // 
            this.button_correct.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.button_correct.Location = new System.Drawing.Point(751, 589);
            this.button_correct.Name = "button_correct";
            this.button_correct.Size = new System.Drawing.Size(95, 23);
            this.button_correct.TabIndex = 8;
            this.button_correct.Text = "Correct!";
            this.button_correct.UseVisualStyleBackColor = true;
            this.button_correct.MouseClick += new System.Windows.Forms.MouseEventHandler(this.button_correct_MouseClick);
            // 
            // groupbox_gridSearchParameters
            // 
            this.groupbox_gridSearchParameters.Anchor = System.Windows.Forms.AnchorStyles.Left;
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_blockSize);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_precision);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_b3Multilplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_b2Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_b1Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_B0Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_a3Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_a2Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_a1Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.nud_A0Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_blockSize);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_precision);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_b3Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_b2Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_b1Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_B0Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_a3Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_a2Multiplier);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_a1Max);
            this.groupbox_gridSearchParameters.Controls.Add(this.label_A0Max);
            this.groupbox_gridSearchParameters.Location = new System.Drawing.Point(12, 383);
            this.groupbox_gridSearchParameters.Name = "groupbox_gridSearchParameters";
            this.groupbox_gridSearchParameters.Size = new System.Drawing.Size(388, 187);
            this.groupbox_gridSearchParameters.TabIndex = 5;
            this.groupbox_gridSearchParameters.TabStop = false;
            this.groupbox_gridSearchParameters.Text = "Grid Search Parameters";
            // 
            // nud_blockSize
            // 
            this.nud_blockSize.Location = new System.Drawing.Point(207, 158);
            this.nud_blockSize.Maximum = new decimal(new int[] {
            20,
            0,
            0,
            0});
            this.nud_blockSize.Name = "nud_blockSize";
            this.nud_blockSize.Size = new System.Drawing.Size(105, 20);
            this.nud_blockSize.TabIndex = 10;
            // 
            // nud_precision
            // 
            this.nud_precision.Location = new System.Drawing.Point(207, 132);
            this.nud_precision.Maximum = new decimal(new int[] {
            8,
            0,
            0,
            0});
            this.nud_precision.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_precision.Name = "nud_precision";
            this.nud_precision.Size = new System.Drawing.Size(105, 20);
            this.nud_precision.TabIndex = 9;
            this.nud_precision.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            // 
            // nud_b3Multilplier
            // 
            this.nud_b3Multilplier.DecimalPlaces = 4;
            this.nud_b3Multilplier.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_b3Multilplier.Location = new System.Drawing.Point(274, 101);
            this.nud_b3Multilplier.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_b3Multilplier.Name = "nud_b3Multilplier";
            this.nud_b3Multilplier.Size = new System.Drawing.Size(105, 20);
            this.nud_b3Multilplier.TabIndex = 8;
            this.nud_b3Multilplier.Value = new decimal(new int[] {
            4,
            0,
            0,
            131072});
            // 
            // nud_b2Multiplier
            // 
            this.nud_b2Multiplier.DecimalPlaces = 4;
            this.nud_b2Multiplier.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_b2Multiplier.Location = new System.Drawing.Point(274, 75);
            this.nud_b2Multiplier.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_b2Multiplier.Name = "nud_b2Multiplier";
            this.nud_b2Multiplier.Size = new System.Drawing.Size(105, 20);
            this.nud_b2Multiplier.TabIndex = 6;
            this.nud_b2Multiplier.Value = new decimal(new int[] {
            4,
            0,
            0,
            131072});
            // 
            // nud_b1Max
            // 
            this.nud_b1Max.DecimalPlaces = 3;
            this.nud_b1Max.Increment = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            this.nud_b1Max.Location = new System.Drawing.Point(274, 49);
            this.nud_b1Max.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_b1Max.Name = "nud_b1Max";
            this.nud_b1Max.Size = new System.Drawing.Size(105, 20);
            this.nud_b1Max.TabIndex = 4;
            this.nud_b1Max.Value = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            // 
            // nud_B0Max
            // 
            this.nud_B0Max.Location = new System.Drawing.Point(274, 23);
            this.nud_B0Max.Maximum = new decimal(new int[] {
            20,
            0,
            0,
            0});
            this.nud_B0Max.Name = "nud_B0Max";
            this.nud_B0Max.Size = new System.Drawing.Size(105, 20);
            this.nud_B0Max.TabIndex = 2;
            this.nud_B0Max.Value = new decimal(new int[] {
            3,
            0,
            0,
            0});
            // 
            // nud_a3Multiplier
            // 
            this.nud_a3Multiplier.DecimalPlaces = 4;
            this.nud_a3Multiplier.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_a3Multiplier.Location = new System.Drawing.Point(86, 101);
            this.nud_a3Multiplier.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_a3Multiplier.Name = "nud_a3Multiplier";
            this.nud_a3Multiplier.Size = new System.Drawing.Size(105, 20);
            this.nud_a3Multiplier.TabIndex = 7;
            this.nud_a3Multiplier.Value = new decimal(new int[] {
            4,
            0,
            0,
            131072});
            // 
            // nud_a2Multiplier
            // 
            this.nud_a2Multiplier.DecimalPlaces = 4;
            this.nud_a2Multiplier.Increment = new decimal(new int[] {
            1,
            0,
            0,
            196608});
            this.nud_a2Multiplier.Location = new System.Drawing.Point(86, 75);
            this.nud_a2Multiplier.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_a2Multiplier.Name = "nud_a2Multiplier";
            this.nud_a2Multiplier.Size = new System.Drawing.Size(105, 20);
            this.nud_a2Multiplier.TabIndex = 5;
            this.nud_a2Multiplier.Value = new decimal(new int[] {
            4,
            0,
            0,
            131072});
            // 
            // nud_a1Max
            // 
            this.nud_a1Max.DecimalPlaces = 3;
            this.nud_a1Max.Increment = new decimal(new int[] {
            5,
            0,
            0,
            131072});
            this.nud_a1Max.Location = new System.Drawing.Point(86, 49);
            this.nud_a1Max.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_a1Max.Name = "nud_a1Max";
            this.nud_a1Max.Size = new System.Drawing.Size(105, 20);
            this.nud_a1Max.TabIndex = 3;
            this.nud_a1Max.Value = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            // 
            // nud_A0Max
            // 
            this.nud_A0Max.Location = new System.Drawing.Point(86, 21);
            this.nud_A0Max.Maximum = new decimal(new int[] {
            20,
            0,
            0,
            0});
            this.nud_A0Max.Name = "nud_A0Max";
            this.nud_A0Max.Size = new System.Drawing.Size(105, 20);
            this.nud_A0Max.TabIndex = 1;
            this.nud_A0Max.Value = new decimal(new int[] {
            3,
            0,
            0,
            0});
            // 
            // label_blockSize
            // 
            this.label_blockSize.AutoSize = true;
            this.label_blockSize.Location = new System.Drawing.Point(134, 160);
            this.label_blockSize.Name = "label_blockSize";
            this.label_blockSize.Size = new System.Drawing.Size(57, 13);
            this.label_blockSize.TabIndex = 0;
            this.label_blockSize.Text = "Block Size";
            // 
            // label_precision
            // 
            this.label_precision.AutoSize = true;
            this.label_precision.Location = new System.Drawing.Point(134, 134);
            this.label_precision.Name = "label_precision";
            this.label_precision.Size = new System.Drawing.Size(50, 13);
            this.label_precision.TabIndex = 0;
            this.label_precision.Text = "Precision";
            // 
            // label_b3Multiplier
            // 
            this.label_b3Multiplier.AutoSize = true;
            this.label_b3Multiplier.Location = new System.Drawing.Point(204, 108);
            this.label_b3Multiplier.Name = "label_b3Multiplier";
            this.label_b3Multiplier.Size = new System.Drawing.Size(63, 13);
            this.label_b3Multiplier.TabIndex = 0;
            this.label_b3Multiplier.Text = "b3 Multiplier";
            // 
            // label_b2Multiplier
            // 
            this.label_b2Multiplier.AutoSize = true;
            this.label_b2Multiplier.Location = new System.Drawing.Point(204, 82);
            this.label_b2Multiplier.Name = "label_b2Multiplier";
            this.label_b2Multiplier.Size = new System.Drawing.Size(63, 13);
            this.label_b2Multiplier.TabIndex = 0;
            this.label_b2Multiplier.Text = "b2 Multiplier";
            // 
            // label_b1Max
            // 
            this.label_b1Max.AutoSize = true;
            this.label_b1Max.Location = new System.Drawing.Point(204, 56);
            this.label_b1Max.Name = "label_b1Max";
            this.label_b1Max.Size = new System.Drawing.Size(42, 13);
            this.label_b1Max.TabIndex = 0;
            this.label_b1Max.Text = "b1 Max";
            // 
            // label_B0Max
            // 
            this.label_B0Max.AutoSize = true;
            this.label_B0Max.Location = new System.Drawing.Point(204, 30);
            this.label_B0Max.Name = "label_B0Max";
            this.label_B0Max.Size = new System.Drawing.Size(43, 13);
            this.label_B0Max.TabIndex = 0;
            this.label_B0Max.Text = "B0 Max";
            // 
            // label_a3Multiplier
            // 
            this.label_a3Multiplier.AutoSize = true;
            this.label_a3Multiplier.Location = new System.Drawing.Point(12, 108);
            this.label_a3Multiplier.Name = "label_a3Multiplier";
            this.label_a3Multiplier.Size = new System.Drawing.Size(63, 13);
            this.label_a3Multiplier.TabIndex = 0;
            this.label_a3Multiplier.Text = "a3 Multiplier";
            // 
            // label_a2Multiplier
            // 
            this.label_a2Multiplier.AutoSize = true;
            this.label_a2Multiplier.Location = new System.Drawing.Point(12, 82);
            this.label_a2Multiplier.Name = "label_a2Multiplier";
            this.label_a2Multiplier.Size = new System.Drawing.Size(63, 13);
            this.label_a2Multiplier.TabIndex = 0;
            this.label_a2Multiplier.Text = "a2 Multiplier";
            // 
            // label_a1Max
            // 
            this.label_a1Max.AutoSize = true;
            this.label_a1Max.Location = new System.Drawing.Point(12, 56);
            this.label_a1Max.Name = "label_a1Max";
            this.label_a1Max.Size = new System.Drawing.Size(42, 13);
            this.label_a1Max.TabIndex = 0;
            this.label_a1Max.Text = "a1 Max";
            // 
            // label_A0Max
            // 
            this.label_A0Max.AutoSize = true;
            this.label_A0Max.Location = new System.Drawing.Point(12, 30);
            this.label_A0Max.Name = "label_A0Max";
            this.label_A0Max.Size = new System.Drawing.Size(43, 13);
            this.label_A0Max.TabIndex = 0;
            this.label_A0Max.Text = "A0 Max";
            // 
            // groupbox_simplexParameters
            // 
            this.groupbox_simplexParameters.Anchor = System.Windows.Forms.AnchorStyles.Left;
            this.groupbox_simplexParameters.Controls.Add(this.nud_iterations);
            this.groupbox_simplexParameters.Controls.Add(this.label_iterations);
            this.groupbox_simplexParameters.Controls.Add(this.nud_reflection);
            this.groupbox_simplexParameters.Controls.Add(this.label_reflection);
            this.groupbox_simplexParameters.Controls.Add(this.nud_contraction);
            this.groupbox_simplexParameters.Controls.Add(this.label_contraction);
            this.groupbox_simplexParameters.Controls.Add(this.nud_growth);
            this.groupbox_simplexParameters.Controls.Add(this.label_growth);
            this.groupbox_simplexParameters.Location = new System.Drawing.Point(434, 196);
            this.groupbox_simplexParameters.Name = "groupbox_simplexParameters";
            this.groupbox_simplexParameters.Size = new System.Drawing.Size(197, 131);
            this.groupbox_simplexParameters.TabIndex = 6;
            this.groupbox_simplexParameters.TabStop = false;
            this.groupbox_simplexParameters.Text = "Simplex Parameters";
            // 
            // nud_iterations
            // 
            this.nud_iterations.Increment = new decimal(new int[] {
            100,
            0,
            0,
            0});
            this.nud_iterations.Location = new System.Drawing.Point(86, 97);
            this.nud_iterations.Maximum = new decimal(new int[] {
            50000,
            0,
            0,
            0});
            this.nud_iterations.Minimum = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            this.nud_iterations.Name = "nud_iterations";
            this.nud_iterations.Size = new System.Drawing.Size(105, 20);
            this.nud_iterations.TabIndex = 4;
            this.nud_iterations.Value = new decimal(new int[] {
            10000,
            0,
            0,
            0});
            // 
            // label_iterations
            // 
            this.label_iterations.AutoSize = true;
            this.label_iterations.Location = new System.Drawing.Point(13, 104);
            this.label_iterations.Name = "label_iterations";
            this.label_iterations.Size = new System.Drawing.Size(50, 13);
            this.label_iterations.TabIndex = 0;
            this.label_iterations.Text = "Iterations";
            // 
            // nud_reflection
            // 
            this.nud_reflection.DecimalPlaces = 3;
            this.nud_reflection.Increment = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            this.nud_reflection.Location = new System.Drawing.Point(86, 71);
            this.nud_reflection.Maximum = new decimal(new int[] {
            20,
            0,
            0,
            0});
            this.nud_reflection.Name = "nud_reflection";
            this.nud_reflection.Size = new System.Drawing.Size(105, 20);
            this.nud_reflection.TabIndex = 3;
            this.nud_reflection.Value = new decimal(new int[] {
            13,
            0,
            0,
            65536});
            // 
            // label_reflection
            // 
            this.label_reflection.AutoSize = true;
            this.label_reflection.Location = new System.Drawing.Point(13, 78);
            this.label_reflection.Name = "label_reflection";
            this.label_reflection.Size = new System.Drawing.Size(55, 13);
            this.label_reflection.TabIndex = 0;
            this.label_reflection.Text = "Reflection";
            // 
            // nud_contraction
            // 
            this.nud_contraction.DecimalPlaces = 3;
            this.nud_contraction.Increment = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            this.nud_contraction.Location = new System.Drawing.Point(86, 45);
            this.nud_contraction.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nud_contraction.Name = "nud_contraction";
            this.nud_contraction.Size = new System.Drawing.Size(105, 20);
            this.nud_contraction.TabIndex = 2;
            this.nud_contraction.Value = new decimal(new int[] {
            7,
            0,
            0,
            65536});
            // 
            // label_contraction
            // 
            this.label_contraction.AutoSize = true;
            this.label_contraction.Location = new System.Drawing.Point(13, 52);
            this.label_contraction.Name = "label_contraction";
            this.label_contraction.Size = new System.Drawing.Size(61, 13);
            this.label_contraction.TabIndex = 0;
            this.label_contraction.Text = "Contraction";
            // 
            // nud_growth
            // 
            this.nud_growth.DecimalPlaces = 3;
            this.nud_growth.Increment = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            this.nud_growth.Location = new System.Drawing.Point(86, 19);
            this.nud_growth.Maximum = new decimal(new int[] {
            2,
            0,
            0,
            0});
            this.nud_growth.Name = "nud_growth";
            this.nud_growth.Size = new System.Drawing.Size(105, 20);
            this.nud_growth.TabIndex = 1;
            this.nud_growth.Value = new decimal(new int[] {
            18,
            0,
            0,
            65536});
            // 
            // label_growth
            // 
            this.label_growth.AutoSize = true;
            this.label_growth.Location = new System.Drawing.Point(13, 26);
            this.label_growth.Name = "label_growth";
            this.label_growth.Size = new System.Drawing.Size(41, 13);
            this.label_growth.TabIndex = 0;
            this.label_growth.Text = "Growth";
            // 
            // button_sliverload
            // 
            this.button_sliverload.Location = new System.Drawing.Point(586, 354);
            this.button_sliverload.Name = "button_sliverload";
            this.button_sliverload.Size = new System.Drawing.Size(90, 23);
            this.button_sliverload.TabIndex = 3;
            this.button_sliverload.Text = "Load Sliver";
            this.button_sliverload.UseVisualStyleBackColor = true;
            this.button_sliverload.Click += new System.EventHandler(this.button_sliverload_Click);
            // 
            // button_cleardata
            // 
            this.button_cleardata.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.button_cleardata.Location = new System.Drawing.Point(759, 383);
            this.button_cleardata.Name = "button_cleardata";
            this.button_cleardata.Size = new System.Drawing.Size(90, 23);
            this.button_cleardata.TabIndex = 7;
            this.button_cleardata.Text = "Clear Data";
            this.button_cleardata.UseVisualStyleBackColor = true;
            this.button_cleardata.Click += new System.EventHandler(this.button_cleardata_Click);
            // 
            // label_correctionOrder
            // 
            this.label_correctionOrder.Anchor = System.Windows.Forms.AnchorStyles.Right;
            this.label_correctionOrder.AutoSize = true;
            this.label_correctionOrder.Location = new System.Drawing.Point(413, 404);
            this.label_correctionOrder.Name = "label_correctionOrder";
            this.label_correctionOrder.Size = new System.Drawing.Size(84, 13);
            this.label_correctionOrder.TabIndex = 0;
            this.label_correctionOrder.Text = "Correction Order";
            // 
            // checkedList_layers
            // 
            this.checkedList_layers.CheckOnClick = true;
            this.checkedList_layers.FormattingEnabled = true;
            this.checkedList_layers.Location = new System.Drawing.Point(415, 422);
            this.checkedList_layers.Name = "checkedList_layers";
            this.checkedList_layers.Size = new System.Drawing.Size(174, 139);
            this.checkedList_layers.TabIndex = 0;
            this.checkedList_layers.TabStop = false;
            this.checkedList_layers.ItemCheck += new System.Windows.Forms.ItemCheckEventHandler(this.checkedList_layers_ItemCheck);
            // 
            // button_rotateSliver
            // 
            this.button_rotateSliver.Anchor = System.Windows.Forms.AnchorStyles.Right;
            this.button_rotateSliver.BackgroundImage = global::ArgoGUI.Properties.Resources.rotate_counter_clockwise_512;
            this.button_rotateSliver.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Stretch;
            this.button_rotateSliver.Enabled = false;
            this.button_rotateSliver.Location = new System.Drawing.Point(653, 138);
            this.button_rotateSliver.Name = "button_rotateSliver";
            this.button_rotateSliver.Size = new System.Drawing.Size(23, 23);
            this.button_rotateSliver.TabIndex = 0;
            this.button_rotateSliver.TabStop = false;
            this.button_rotateSliver.UseVisualStyleBackColor = true;
            this.button_rotateSliver.Click += new System.EventHandler(this.button_rotateSliver_Click);
            // 
            // button_flipSliver
            // 
            this.button_flipSliver.Anchor = System.Windows.Forms.AnchorStyles.Right;
            this.button_flipSliver.BackgroundImage = global::ArgoGUI.Properties.Resources.flip_horizontal_object_icone_4623_128;
            this.button_flipSliver.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Stretch;
            this.button_flipSliver.Enabled = false;
            this.button_flipSliver.Location = new System.Drawing.Point(653, 167);
            this.button_flipSliver.Name = "button_flipSliver";
            this.button_flipSliver.Size = new System.Drawing.Size(23, 23);
            this.button_flipSliver.TabIndex = 0;
            this.button_flipSliver.TabStop = false;
            this.button_flipSliver.UseVisualStyleBackColor = true;
            this.button_flipSliver.Click += new System.EventHandler(this.button_flipSliver_Click);
            // 
            // button_rotateBase
            // 
            this.button_rotateBase.Anchor = System.Windows.Forms.AnchorStyles.Right;
            this.button_rotateBase.BackgroundImage = global::ArgoGUI.Properties.Resources.rotate_counter_clockwise_512;
            this.button_rotateBase.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Stretch;
            this.button_rotateBase.Enabled = false;
            this.button_rotateBase.Location = new System.Drawing.Point(428, 138);
            this.button_rotateBase.Name = "button_rotateBase";
            this.button_rotateBase.Size = new System.Drawing.Size(23, 23);
            this.button_rotateBase.TabIndex = 0;
            this.button_rotateBase.TabStop = false;
            this.button_rotateBase.UseVisualStyleBackColor = true;
            this.button_rotateBase.Click += new System.EventHandler(this.button_rotateBase_Click);
            // 
            // button_flipBase
            // 
            this.button_flipBase.Anchor = System.Windows.Forms.AnchorStyles.Right;
            this.button_flipBase.BackgroundImage = global::ArgoGUI.Properties.Resources.flip_horizontal_object_icone_4623_128;
            this.button_flipBase.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Stretch;
            this.button_flipBase.Enabled = false;
            this.button_flipBase.Location = new System.Drawing.Point(428, 167);
            this.button_flipBase.Name = "button_flipBase";
            this.button_flipBase.Size = new System.Drawing.Size(23, 23);
            this.button_flipBase.TabIndex = 0;
            this.button_flipBase.TabStop = false;
            this.button_flipBase.UseVisualStyleBackColor = true;
            this.button_flipBase.Click += new System.EventHandler(this.button_flipBase_Click);
            // 
            // image_sliver
            // 
            this.image_sliver.Location = new System.Drawing.Point(682, 25);
            this.image_sliver.Name = "image_sliver";
            this.image_sliver.Size = new System.Drawing.Size(167, 352);
            this.image_sliver.SizeMode = System.Windows.Forms.PictureBoxSizeMode.CenterImage;
            this.image_sliver.TabIndex = 1;
            this.image_sliver.TabStop = false;
            // 
            // image_base
            // 
            this.image_base.Location = new System.Drawing.Point(12, 25);
            this.image_base.Name = "image_base";
            this.image_base.Size = new System.Drawing.Size(410, 352);
            this.image_base.SizeMode = System.Windows.Forms.PictureBoxSizeMode.CenterImage;
            this.image_base.TabIndex = 0;
            this.image_base.TabStop = false;
            // 
            // button_defaultParameters
            // 
            this.button_defaultParameters.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.button_defaultParameters.Location = new System.Drawing.Point(606, 383);
            this.button_defaultParameters.Name = "button_defaultParameters";
            this.button_defaultParameters.Size = new System.Drawing.Size(123, 23);
            this.button_defaultParameters.TabIndex = 10;
            this.button_defaultParameters.Text = "Default Parameters";
            this.button_defaultParameters.UseVisualStyleBackColor = true;
            this.button_defaultParameters.Click += new System.EventHandler(this.button_defaultParameters_Click);
            // 
            // Panel
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(96F, 96F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Dpi;
            this.AutoSize = true;
            this.ClientSize = new System.Drawing.Size(872, 650);
            this.Controls.Add(this.button_defaultParameters);
            this.Controls.Add(this.button_rotateSliver);
            this.Controls.Add(this.button_flipSliver);
            this.Controls.Add(this.button_rotateBase);
            this.Controls.Add(this.button_flipBase);
            this.Controls.Add(this.checkedList_layers);
            this.Controls.Add(this.label_correctionOrder);
            this.Controls.Add(this.button_cleardata);
            this.Controls.Add(this.button_sliverload);
            this.Controls.Add(this.groupbox_simplexParameters);
            this.Controls.Add(this.groupbox_gridSearchParameters);
            this.Controls.Add(label_sliverImage);
            this.Controls.Add(this.button_correct);
            this.Controls.Add(this.button_baseload);
            this.Controls.Add(this.progressbar);
            this.Controls.Add(this.groupbox_correctionMode);
            this.Controls.Add(label_baseImage);
            this.Controls.Add(this.dropdown_imageLayer);
            this.Controls.Add(label_imageLayer);
            this.Controls.Add(this.image_sliver);
            this.Controls.Add(this.image_base);
            this.Name = "Panel";
            this.StartPosition = System.Windows.Forms.FormStartPosition.Manual;
            this.Text = "AroGUI";
            this.Load += new System.EventHandler(this.Panel_Load);
            this.groupbox_correctionMode.ResumeLayout(false);
            this.groupbox_correctionMode.PerformLayout();
            this.groupbox_gridSearchParameters.ResumeLayout(false);
            this.groupbox_gridSearchParameters.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_blockSize)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_precision)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_b3Multilplier)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_b2Multiplier)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_b1Max)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_B0Max)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_a3Multiplier)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_a2Multiplier)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_a1Max)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_A0Max)).EndInit();
            this.groupbox_simplexParameters.ResumeLayout(false);
            this.groupbox_simplexParameters.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nud_iterations)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_reflection)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_contraction)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nud_growth)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.image_sliver)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.image_base)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }
    }
}

