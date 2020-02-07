using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace pi2cs
{
    /// <summary>
    /// Shows Pi2 image with zoom, scroll etc. controls.
    /// </summary>
    public partial class Pi2PictureViewer : UserControl
    {
        
        /// <summary>
        /// Use this to 
        /// </summary>
        internal Pi2PictureBox PictureBox
        {
            get;
            private set;
        }

        public Pi2PictureToolStrip ToolStrip
        {
            get;
            set;
        }

        /// <summary>
        /// Gets and sets the zoom level.
        /// </summary>
        public float Zoom
        {
            get
            {
                return PictureBox.Zoom;
            }
            set
            {
                PictureBox.Zoom = value;
                UpdateUI();
            }
        }

        /// <summary>
        /// Gray value that is mapped to black.
        /// </summary>
        public float Min
        {
            get
            {
                return PictureBox.Min;
            }
            set
            {
                PictureBox.Min = value;
                UpdateUI();
            }
        }

        /// <summary>
        /// Gray value that is mapped to white.
        /// </summary>
        public float Max
        {
            get
            {
                return PictureBox.Max;
            }
            set
            {
                PictureBox.Max = value;
                UpdateUI();
            }
        }

        /// <summary>
        /// Minimum pixel value in the original image.
        /// </summary>
        public float GlobalMin
        {
            get
            {
                return PictureBox.GlobalMin;
            }
        }

        /// <summary>
        /// Maximum pixel value in the original image.
        /// </summary>
        public float GlobalMax
        {
            get
            {
                return PictureBox.GlobalMax;
            }
        }

        /// <summary>
        /// Minimum possible value in the original image.
        /// </summary>
        public float DynamicMin
        {
            get
            {
                return PictureBox.DynamicMin;
            }
        }

        /// <summary>
        /// Maximum possible value in the original image.
        /// </summary>
        public float DynamicMax
        {
            get
            {
                return PictureBox.DynamicMax;
            }
        }

        /// <summary>
        /// Pi image to show.
        /// </summary>
        public Pi2Image PiImage
        {
            get
            {
                return PictureBox.PiImage;
            }
            set
            {
                PictureBox.PiImage = value;
                PictureBox.UpdateImage();
                //SetAutoContrast();
            }
        }

        
        /// <summary>
        /// Profile window.
        /// </summary>
        private ProfileDialog Profile;

        /// <summary>
        /// Constructor
        /// </summary>
        public Pi2PictureViewer()
        {
            InitializeComponent();
            
            labelPixelValue.Text = String.Empty;
            labelSliceNumber.Text = String.Empty;
            scrollBarSlice.Visible = false;

            PictureBox = new Pi2PictureBox();
            PictureBox.AutoResize = false;
            PictureBox.Dock = DockStyle.Fill;
            PictureBox.MouseDown += PictureBox_MouseDown;
            PictureBox.MouseMove += PictureBox_MouseMove;
            PictureBox.MouseUp += PictureBox_MouseUp;
            PictureBox.MouseWheel += PictureBox_MouseWheel;
            Controls.Add(PictureBox);

            //Mode = MouseMode.AddAnnotation;
            Mode = MouseMode.Pan;

            Profile = new ProfileDialog(this);


            // DEBUG: Test annotations
            //NewAnnotationName = "Line segment";
            NewAnnotationName = "Horizontal line";

            //HLineAnnotation ann = new HLineAnnotation();
            //ann.ControlPoints.Add(new Vec3(0, 60, 0));
            //PictureBox.Annotations.Add(ann);

            //LineSegmentAnnotation ann2 = new LineSegmentAnnotation();
            //ann2.ControlPoints.Add(new Vec3(10, 20, 0));
            //ann2.ControlPoints.Add(new Vec3(200, 180, 0));
            //ann2.LineWidth = 2;
            //ann2.DrawToAllSlices = true;
            //ann2.Color = Color.Red;
            //PictureBox.Annotations.Add(ann2);
        }
        
        /// <summary>
        /// Sets contrast of the image automatically.
        /// </summary>
        public void SetAutoContrast()
        {
            Min = PictureBox.GlobalMin;
            Max = PictureBox.GlobalMax;
        }

        /// <summary>
        /// Zooms so that the image fits into the window.
        /// Requires tool strip.
        /// </summary>
        public void ZoomToFit()
        {
            if (ToolStrip != null)
                ToolStrip.ZoomFit();
        }


        #region Profile plot

        internal Annotation ProfileAnnotation;

        /// <summary>
        /// Removes profile annotation from the image.
        /// </summary>
        internal void EraseProfileAnnotation()
        {
            if (ProfileAnnotation != null)
            {
                PictureBox.Annotations.Remove(ProfileAnnotation);
                ProfileAnnotation = null;
                Refresh();
            }
        }

        /// <summary>
        /// Updates profile plot if profile annotation has been made.
        /// </summary>
        internal void UpdateProfile()
        {
            if (Mode == MouseMode.Profile && ProfileAnnotation != null && ProfileAnnotation.ControlPoints.Count >= 2)
            {
                Vec2 profileStart = ProfileAnnotation.ControlPoints[0].GetXY();
                Vec2 profileEnd = ProfileAnnotation.ControlPoints[1].GetXY();
                ProfileAnnotation.LineWidth = Profile.ProfileWidth;
                double profileRadius = ProfileAnnotation.LineWidth / 2;
                
                Vec2 delta = profileEnd - profileStart;
                double l = delta.Norm();
                Vec2 n = new Vec2(delta.Y, -delta.X) / l; // Unit normal

                // Extract profile data
                Profile.Clear();
                double tstep = 1.0 / Math.Max(Math.Abs(delta.X), Math.Abs(delta.Y));
                for (double t = 0; t <= 1; t += tstep)
                {
                    Vec2 f = profileStart + t * delta;
                    double r = t * l;

                    double sum = 0.0;
                    double count = 0.0;

                    for (double s = 0; s < profileRadius; s++)
                    {
                        sum += PictureBox.GetValueImage(f + s * n);
                        count++;

                        if (s > 0)
                        {
                            sum += PictureBox.GetValueImage(f - s * n);
                            count++;
                        }
                    }

                    Profile.AddPoint(r, sum / count);
                }

                Profile.Location = PointToScreen(Point.Round(PictureBox.PictureToScreen(profileEnd).ToPointF()));
                Profile.Show();
            }
        }
        
        #endregion

        

        internal enum MouseMode
        {
            /// <summary>
            /// Mouse picks existing annotations and control points.
            /// </summary>
            Pick,
            /// <summary>
            /// Mouse pans the image.
            /// </summary>
            Pan,
            /// <summary>
            /// Mouse adds annotations.
            /// </summary>
            AddAnnotation,
            /// <summary>
            /// Delete clicked annotations
            /// </summary>
            Erase,
            /// <summary>
            /// Plot profile
            /// </summary>
            Profile
        }

        private MouseMode mode;

        /// <summary>
        /// Determines how the user can interact with the image using mouse cursor.
        /// </summary>
        internal MouseMode Mode
        {
            get
            {
                return mode;
            }
            set
            {
                mode = value;

                if (mode == MouseMode.Pan)
                    PictureBox.Cursor = Cursors.Hand;
                else if (mode == MouseMode.Pick || mode == MouseMode.Erase)
                    PictureBox.Cursor = Cursors.Arrow;
                else
                    PictureBox.Cursor = Cursors.Cross;
            }
        }

        /// <summary>
        /// Adds annotation that requires no control points.
        /// There can be only one of those annotations per annotation name.
        /// </summary>
        /// <param name="annotationName"></param>
        public void AddSimpleAnnotation(string annotationName)
        {
            Annotation ann = AnnotationFactory.Create(annotationName);

            foreach(Annotation ann2 in PictureBox.Annotations)
            {
                if (ann2.GetType() == ann.GetType())
                    return; // This annotation is already visible, do not add a second one.
            }

            PictureBox.Annotations.Add(ann);
            Refresh();
        }

        /// <summary>
        /// Mouse position when pan operation started.
        /// </summary>
        private Point PanStartMousePosition = Point.Empty;

        /// <summary>
        /// Picture position when pan operation started.
        /// </summary>
        private PointF PanStartPicturePosition;

        /// <summary>
        /// Annotation that is currently active.
        /// </summary>
        private Annotation ActiveAnnotation;

        /// <summary>
        /// Index of active control point or -1.
        /// </summary>
        private int ActiveControlPoint;

        /// <summary>
        /// Name of next annotation to add.
        /// </summary>
        public string NewAnnotationName;

        /// <summary>
        /// Stores old locations of control points while moving them.
        /// </summary>
        public List<Vec3> oldControlPointLocations = new List<Vec3>();

        /// <summary>
        /// Mouse position where control point move started.
        /// </summary>
        public Vec3 controlPointSelectLocation;


        private void PictureBox_MouseDown(object sender, MouseEventArgs e)
        {
            if((Mode == MouseMode.Pan && e.Button == MouseButtons.Left) ||
                (e.Button == MouseButtons.Middle))
            {
                // Start panning the image
                PanStartMousePosition = e.Location;
                PanStartPicturePosition = PictureBox.PicturePosition;
            }
            else if((Mode == MouseMode.Pick || Mode == MouseMode.Erase) &&
                e.Button == MouseButtons.Left &&
                ActiveAnnotation == null)
            {
                // Pick and activate annotation and possibly also control point

                var result = PictureBox.FindAnnotation(e.Location);

                if(result != null)
                {
                    if (Mode == MouseMode.Pick)
                    {
                        ActiveAnnotation = result.Item1;
                        ActiveControlPoint = result.Item2;

                        oldControlPointLocations.Clear();
                        oldControlPointLocations.AddRange(ActiveAnnotation.ControlPoints);

                        controlPointSelectLocation = new Vec3(PictureBox.ScreenToPicture(new Vec2(e.Location)), PictureBox.Slice);
                    }
                    else
                    {
                        // Erase annotation
                        PictureBox.Annotations.Remove(result.Item1);
                        Refresh();
                    }
                }
            }
            else if((Mode == MouseMode.AddAnnotation || Mode == MouseMode.Profile) && 
                e.Button == MouseButtons.Left &&
                ActiveAnnotation == null)
            {
                if (Mode == MouseMode.Profile)
                {
                    NewAnnotationName = "Line segment";
                    EraseProfileAnnotation();
                }

                // Add new annotation object
                if(!String.IsNullOrEmpty(NewAnnotationName))
                {
                    ActiveAnnotation = AnnotationFactory.Create(NewAnnotationName);
                    PictureBox.Annotations.Add(ActiveAnnotation);
                    Vec3 pos = new Vec3(PictureBox.ScreenToPicture(new Vec2(e.Location)), PictureBox.Slice);
                    ActiveAnnotation.ControlPoints.Add(pos);
                    ActiveControlPoint = 0;

                    // Add another control point right away if multiple points are required.
                    if(ActiveAnnotation.RequiredPointCount > 1)
                    {
                        ActiveAnnotation.ControlPoints.Add(pos);
                        ActiveControlPoint = 1;
                    }

                    Refresh();

                    if(Mode != MouseMode.Profile)
                        if(ToolStrip != null)
                            ToolStrip.ResetAnnotation();
                }

                if (Mode == MouseMode.Profile)
                {
                    ProfileAnnotation = ActiveAnnotation;
                    UpdateProfile();
                }
            }
        }
        
        private void PictureBox_MouseMove(object sender, MouseEventArgs e)
        {
            // Update pixel value label
            float val = PictureBox.GetValue(new Vec2(e.Location));
            string valstr;
            if (float.IsNaN(val))
                valstr = "No value";
            else
                valstr = val.ToString();

            Vec2 picPos = PictureBox.ScreenToPicture(new Vec2(e.Location));
            if(PictureBox.Zoom != 1)
                labelPixelValue.Text = $"({picPos.X:F2}, {picPos.Y:F2}): {valstr}";
            else
                labelPixelValue.Text = $"({picPos.X:F0}, {picPos.Y:F0}): {valstr}";
            
            if(PanStartMousePosition != Point.Empty)
            {
                // Pan image

                Vec2 delta = (new Vec2(e.Location) - new Vec2(PanStartMousePosition)) / Zoom;
                PictureBox.PicturePosition = (new Vec2(PanStartPicturePosition) - delta).ToPointF();

                UpdateUI();
            }
            else if(ActiveAnnotation != null && e.Button == MouseButtons.Left)
            {
                if(ActiveControlPoint >= 0)
                {
                    // Move control point
                    Vec3 pos = new Vec3(PictureBox.ScreenToPicture(new Vec2(e.Location)), PictureBox.Slice);
                    ActiveAnnotation.ControlPoints[ActiveControlPoint] = pos;
                    
                }
                else
                {
                    // Move all control points
                    Vec3 delta = new Vec3(PictureBox.ScreenToPicture(new Vec2(e.Location)), PictureBox.Slice) - controlPointSelectLocation;

                    for (int n = 0; n < ActiveAnnotation.ControlPoints.Count; n++)
                        ActiveAnnotation.ControlPoints[n] = oldControlPointLocations[n] + delta;
                }
                Refresh();

                UpdateProfile();
            }
        }

        private void PictureBox_MouseUp(object sender, MouseEventArgs e)
        {
            if ((Mode == MouseMode.Pan && e.Button == MouseButtons.Left) ||
                (e.Button == MouseButtons.Middle))
            {
                // Stop panning
                PanStartMousePosition = Point.Empty;
            }
            else if(ActiveAnnotation != null && ActiveControlPoint >= 0 && e.Button == MouseButtons.Left)
            {
                // Move current control point to the final location of button up event.
                Vec3 pos = new Vec3(PictureBox.ScreenToPicture(new Vec2(e.Location)), PictureBox.Slice);
                ActiveAnnotation.ControlPoints[ActiveControlPoint] = pos;

                if (ActiveAnnotation.ControlPoints.Count < ActiveAnnotation.RequiredPointCount)
                {
                    // Add next control point
                    ActiveAnnotation.ControlPoints.Add(pos);
                    ActiveControlPoint = ActiveAnnotation.ControlPoints.Count - 1;
                    Refresh();
                }
                else
                {
                    // Disable control point move mode.
                    ActiveAnnotation = null;
                    ActiveControlPoint = -1;
                }

                UpdateProfile();
            }
            else if(ActiveAnnotation != null && e.Button == MouseButtons.Left)
            {
                // End move of all control points
                ActiveAnnotation = null;
            }
        }

        /// <summary>
        /// Scrolls slices
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void PictureBox_MouseWheel(object sender, MouseEventArgs e)
        {
            int delta = -Math.Sign(e.Delta);

            if ((Control.ModifierKeys & Keys.Control) != 0)
            {
                // Zoom
                if(ToolStrip != null)
                {
                    if (delta > 0)
                        ToolStrip.ZoomOut();
                    else
                        ToolStrip.ZoomIn();
                }
            }
            else
            {
                // Change slice
                int newVal = scrollBarSlice.Value + delta;
                if (newVal < scrollBarSlice.Minimum)
                    newVal = scrollBarSlice.Minimum;
                if (newVal > scrollBarSlice.Maximum - scrollBarSlice.LargeChange + 1)
                    newVal = scrollBarSlice.Maximum - scrollBarSlice.LargeChange + 1;
                scrollBarSlice.Value = newVal;
            }
        }

        /// <summary>
        /// Updates the user interface.
        /// </summary>
        public void UpdateUI(bool updateImageFromPi = false)
        {
            if (updateImageFromPi)
                PictureBox.UpdateImage();
            else
                PictureBox.UpdateScreen();

            scrollBarSlice.Visible = PictureBox.OriginalDepth > 1;
            scrollBarSlice.Maximum = PictureBox.OriginalDepth + scrollBarSlice.LargeChange - 2; // Note: -2 makes the scroll bar to stop to correct position - probably.
            labelSliceNumber.Text = $"{PictureBox.Slice + 1} / {PictureBox.OriginalDepth}";
        }
        
        private void scrollBarSlice_ValueChanged(object sender, EventArgs e)
        {
            PictureBox.Slice = scrollBarSlice.Value;
            UpdateUI(true);
            UpdateProfile();
        }
    }
}
