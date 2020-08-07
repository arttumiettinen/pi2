using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using pi2cs.Properties;

namespace pi2cs
{
    /// <summary>
    /// ToolStrip control that controls Pi2PictureViewer.
    /// </summary>
    public partial class Pi2PictureToolStrip : ToolStrip
    {
        /// <summary>
        /// String that should be shown in about box in any application using this control.
        /// </summary>
        public static string About
        {
            get
            {
                return "Zoom in and zoom out icons created by icons8 and distributed unders Creative Commons Attribution-NoDerivs 3.0 license (https://www.iconsdb.com/black-icons/zoom-in-icon.html, https://www.iconsdb.com/black-icons/zoom-out-icon.html). " +
                    "Histogram icon is in public domain (https://www.iconsdb.com/black-icons/bar-chart-5-icon.html). " +
                    "Auto-contrast icon is provided by Iconic and licensed under the MIT license (https://www.iconsdb.com/black-icons/contrast-icon.html). " +
                    "Pan icon is provided by icons8 under the Creative Commons Attribution-NoDerivs 3.0 license (https://www.iconsdb.com/black-icons/so-so-icon.html). " +
                    "Eraser icon is provided by icons8 under the Creative Commons Attribution-NoDerivs 3.0 license (https://www.iconsdb.com/black-icons/erase-icon.html). " +
                    "Cursor icon is provided by icons8 under the Creative Commons Attribution-NoDerivs 3.0 license (https://www.iconsdb.com/black-icons/cursor-icon.html). " +
                    "Profile plot icon is provided by icons8 under the Creative Commons Attribution-NoDerivs 3.0 license (https://www.iconsdb.com/black-icons/scatter-plot-icon.html). " +
                    "Line segment icon is provided by icons8 under the Creative Commons Attribution-NoDerivs 3.0 license (https://www.iconsdb.com/black-icons/line-2-icon.html). " +
                    "Sticky annotations icon is provided by Coskun Deniz as CC0 1.0 Universal (CC0 1.0) Public Domain Dedication (https://www.iconsdb.com/black-icons/magnet-icon.html). " +
                    "Zoom to fit icon is provided by Icons8 and provided as linkware (http://icons8.com).";

            }
        }

        
        private Pi2PictureViewer pictureViewer;

        /// <summary>
        /// Gets and sets picture viewer associated with this tool strip.
        /// </summary>
        public Pi2PictureViewer PictureViewer
        {
            get
            {
                return pictureViewer;
            }
            set
            {
                if(pictureViewer != null)
                    pictureViewer.ToolStrip = null;

                pictureViewer = value;

                if(pictureViewer != null)
                    pictureViewer.ToolStrip = this;

                UpdateToolButtons();
            }
        }

        /// <summary>
        /// Histogram window.
        /// </summary>
        private HistogramDialog Histogram;

        /// <summary>
        /// Button that shows histogram.
        /// </summary>
        private ToolStripButton buttonHistogram;

        /// <summary>
        /// Label that shows current zoom level.
        /// </summary>
        private ToolStripLabel labelZoom;

        /// <summary>
        /// Annotation and mouse mode buttons.
        /// </summary>
        private ToolStripButton buttonCross;
        private ToolStripButton buttonErase;
        private ToolStripButton buttonHand;
        private ToolStripButton buttonProfile;
        private List<ToolStripButton> annotationButtons;
        private ToolStripButton buttonStickyAnnotationTools;

        /// <summary>
        /// Constructor
        /// </summary>
        public Pi2PictureToolStrip()
        {
            InitializeComponent();

            Histogram = new HistogramDialog(null);

            ImageScalingSize = new Size(20, 20);

            ToolStripButton buttonZoomOut = new ToolStripButton(String.Empty, Resources.zoom_out_multi_size.ToBitmap(), buttonZoomOutClick);
            buttonZoomOut.ToolTipText = "Zoom out" + Environment.NewLine + "Alternatively use Control + mouse wheel to zoom.";
            Items.Add(buttonZoomOut);

            labelZoom = new ToolStripLabel("100 %");
            labelZoom.ToolTipText = "Current zoom level." + Environment.NewLine + "Use Control + mouse wheel or toolbar buttons to zoom."; ;
            Items.Add(labelZoom);

            ToolStripButton buttonZoomIn = new ToolStripButton(String.Empty, Resources.zoom_in_multi_size.ToBitmap(), buttonZoomInClick);
            buttonZoomIn.ToolTipText = "Zoom in" + Environment.NewLine + "Alternatively use Control + mouse wheel to zoom.";
            Items.Add(buttonZoomIn);

            ToolStripButton buttonZoomFit = new ToolStripButton(String.Empty, Resources.zoom_to_fit.ToBitmap(), buttonZoomFitClick);
            buttonZoomFit.ToolTipText = "Zoom such that the whole image fits into the window.";
            Items.Add(buttonZoomFit);

            ToolStripButton buttonZoomFull = new ToolStripButton(String.Empty, Resources.zoom_100.ToBitmap(), buttonZoomFullClick);
            buttonZoomFull.ToolTipText = "Zooms such that one screen pixel equals one image pixel.";
            Items.Add(buttonZoomFull);


            Items.Add(new ToolStripSeparator());

            buttonHistogram = new ToolStripButton(String.Empty, Resources.bar_chart_5_multi_size.ToBitmap(), buttonHistogramClick);
            buttonHistogram.ToolTipText = "Show histogram";
            Items.Add(buttonHistogram);

            ToolStripButton buttonAutoContrast = new ToolStripButton(String.Empty, Resources.contrast_multi_size.ToBitmap(), buttonAutoContrastClick);
            buttonAutoContrast.ToolTipText = "Adjust contrast automatically";
            Items.Add(buttonAutoContrast);

            Items.Add(new ToolStripSeparator());

            buttonCross = new ToolStripButton(String.Empty, Resources.cursor_multi_size.ToBitmap(), buttonCrossClick);
            buttonCross.ToolTipText = "Move and adjust annotations.";
            Items.Add(buttonCross);

            buttonErase = new ToolStripButton(String.Empty, Resources.erase_multi_size.ToBitmap(), buttonEraseClick);
            buttonErase.ToolTipText = "Erase annotations.";
            Items.Add(buttonErase);

            buttonHand = new ToolStripButton(String.Empty, Resources.so_so_multi_size.ToBitmap(), buttonHandClick);
            buttonHand.ToolTipText = "Pan the image." + Environment.NewLine + "Alternatively, use the middle mouse button to pan the image.";
            Items.Add(buttonHand);

            buttonProfile = new ToolStripButton(String.Empty, Resources.scatter_plot_multi_size.ToBitmap(), buttonProfileClick);
            buttonProfile.ToolTipText = "Plot a profile along a line.";
            Items.Add(buttonProfile);

            Items.Add(new ToolStripSeparator());

            buttonStickyAnnotationTools = new ToolStripButton(String.Empty, Resources.magnet_multi_size.ToBitmap());
            buttonStickyAnnotationTools.ToolTipText = "Activate this button to keep selected annotation tool active indefinitely." + Environment.NewLine + "Deactivate this button to automatically switch to Move and adjust tool after creation of an annotation.";
            buttonStickyAnnotationTools.CheckOnClick = true;
            Items.Add(buttonStickyAnnotationTools);

            Items.Add(new ToolStripSeparator());

            annotationButtons = new List<ToolStripButton>();
            List<string> annotations = AnnotationFactory.Names();
            foreach(string name in annotations)
            {
                ToolStripButton btn = new ToolStripButton(String.Empty, AnnotationFactory.GetIcon(name), buttonActivateAnnotationClick, name);
                btn.ToolTipText = "Create new " + name.ToLower() + " annotation.";
                Items.Add(btn);
                annotationButtons.Add(btn);
            }

        }

        private void buttonActivateAnnotationClick(object sender, EventArgs e)
        {
            if (PictureViewer == null)
                return;

            string annotationName = ((ToolStripButton)sender).Name;

            if (AnnotationFactory.Create(annotationName).RequiredPointCount > 0)
            {
                PictureViewer.NewAnnotationName = annotationName;
                PictureViewer.Mode = Pi2PictureViewer.MouseMode.AddAnnotation;
            }
            else
            {
                PictureViewer.AddSimpleAnnotation(annotationName);
            }
            
            UpdateToolButtons();
        }

        private void buttonCrossClick(object sender, EventArgs e)
        {
            if (PictureViewer == null)
                return;

            PictureViewer.Mode = Pi2PictureViewer.MouseMode.Pick;
            UpdateToolButtons();
        }

        private void buttonEraseClick(object sender, EventArgs e)
        {
            if (PictureViewer == null)
                return;

            PictureViewer.Mode = Pi2PictureViewer.MouseMode.Erase;
            UpdateToolButtons();
        }

        private void buttonHandClick(object sender, EventArgs e)
        {
            if (PictureViewer == null)
                return;

            PictureViewer.Mode = Pi2PictureViewer.MouseMode.Pan;
            UpdateToolButtons();
        }

        private void buttonProfileClick(object sender, EventArgs e)
        {
            if (PictureViewer == null)
                return;

            PictureViewer.Mode = Pi2PictureViewer.MouseMode.Profile;
            UpdateToolButtons();
        }

        /// <summary>
        /// Updates pressed/up state of annotation and mouse mode buttons and other controls on the tool strip.
        /// </summary>
        public void UpdateToolButtons()
        {
            if (PictureViewer != null)
            {
                buttonCross.Checked = PictureViewer.Mode == Pi2PictureViewer.MouseMode.Pick;
                buttonErase.Checked = PictureViewer.Mode == Pi2PictureViewer.MouseMode.Erase;
                buttonHand.Checked = PictureViewer.Mode == Pi2PictureViewer.MouseMode.Pan;
                buttonProfile.Checked = PictureViewer.Mode == Pi2PictureViewer.MouseMode.Profile;

                foreach (ToolStripButton btn in annotationButtons)
                    btn.Checked = PictureViewer.Mode == Pi2PictureViewer.MouseMode.AddAnnotation && String.Equals(btn.Name, PictureViewer.NewAnnotationName);
            }

            UpdateZoomLabel();
        }

        /// <summary>
        /// Sets tool to annotation mover.
        /// </summary>
        public void ResetAnnotation()
        {
            if (!buttonStickyAnnotationTools.Checked)
                buttonCross.PerformClick();
        }

        /// <summary>
        /// Updates text shown in zoom level label.
        /// </summary>
        private void UpdateZoomLabel()
        {
            float zoom = 1;
            if (PictureViewer != null)
                zoom = PictureViewer.Zoom;

            labelZoom.Text = $"{(zoom * 100):F0} %";
        }


        private void buttonZoomFitClick(object sender, EventArgs e)
        {
            if(PictureViewer != null)
                PictureViewer.ZoomFit();
        }

        private void buttonZoomFullClick(object sender, EventArgs e)
        {
            if (PictureViewer != null)
                PictureViewer.ZoomFull();
        }

        private void buttonZoomOutClick(object sender, EventArgs e)
        {
            if (PictureViewer != null)
                PictureViewer.ZoomOut();
        }

        private void buttonZoomInClick(object sender, EventArgs e)
        {
            if (PictureViewer != null)
                PictureViewer.ZoomIn();
        }

        private void buttonHistogramClick(object sender, EventArgs e)
        {
            if (PictureViewer != null)
            {
                Histogram.PictureBox = PictureViewer.PictureBox;
                Histogram.Location = PointToScreen(buttonHistogram.Bounds.Location);
                Histogram.Show();
            }
        }

        private void buttonAutoContrastClick(object sender, EventArgs e)
        {
            if(PictureViewer != null)
                PictureViewer.SetAutoContrast();
        }

    }
}
