using pi2cs.Properties;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
    interface IAnnotationDefinition
    {
        Image Icon
        {
            get;
        }

        Annotation Create();
    }

    class AnnotationDefinition<T> : IAnnotationDefinition where T : Annotation, new()
    {

        public Image Icon
        {
            get;
            private set;
        }

        public Annotation Create()
        {
            return new T();
        }

        public AnnotationDefinition(Image icon)
        {
            Icon = icon;
        }
    }

    /// <summary>
    /// Creates annotation objects.
    /// </summary>
    static class AnnotationFactory
    {
        private static Dictionary<string, IAnnotationDefinition> annotationTypes;

        static AnnotationFactory()
        {
            annotationTypes = new Dictionary<string, IAnnotationDefinition>();
            annotationTypes.Add("Horizontal line", new AnnotationDefinition<HLineAnnotation>(Resources.hline.ToBitmap()));
            Image vline = Resources.hline.ToBitmap();
            vline.RotateFlip(RotateFlipType.Rotate90FlipNone);
            annotationTypes.Add("Vertical line", new AnnotationDefinition<VLineAnnotation>(vline));
            annotationTypes.Add("Line segment", new AnnotationDefinition<LineSegmentAnnotation>(Resources.line_2_multi_size.ToBitmap()));
            annotationTypes.Add("Center cross", new AnnotationDefinition<CenterCrossAnnotation>(Resources.cross.ToBitmap()));
            annotationTypes.Add("Rectangle", new AnnotationDefinition<RectangleAnnotation>(Resources.rectangle.ToBitmap()));
        }

        /// <summary>
        /// Creates annotation, given its name.
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public static Annotation Create(string name)
        {
            return annotationTypes[name].Create();
        }

        /// <summary>
        /// Gets icon for given annotation.
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public static Image GetIcon(string name)
        {
            return annotationTypes[name].Icon;
        }

        /// <summary>
        /// Get a list of names of all annotations.
        /// </summary>
        /// <returns></returns>
        public static List<string> Names()
        {
            return annotationTypes.Keys.ToList();
        }
    }
}
