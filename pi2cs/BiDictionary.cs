using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace pi2cs
{
    /// <summary>
    /// 2-directional dictionary.
    /// </summary>
    /// <typeparam name="T1"></typeparam>
    /// <typeparam name="T2"></typeparam>
    class BiDictionary<T1, T2>
    {
        private Dictionary<T1, T2> forward = new Dictionary<T1, T2>();
        private Dictionary<T2, T1> reverse = new Dictionary<T2, T1>();

        /// <summary>
        /// Adds element to the map.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        public void Add(T1 a, T2 b)
        {
            forward.Add(a, b);
            reverse.Add(b, a);
        }

        /// <summary>
        /// Clears the map.
        /// </summary>
        public void Clear()
        {
            forward.Clear();
            reverse.Clear();
        }

        /// <summary>
        /// Tests if the given key is in the map.
        /// </summary>
        /// <param name="key"></param>
        /// <returns></returns>
        public bool Contains(T1 key)
        {
            return forward.ContainsKey(key);
        }

        /// <summary>
        /// Tests if the given key is in the map.
        /// </summary>
        /// <param name="key"></param>
        /// <returns></returns>
        public bool Contains(T2 key)
        {
            return reverse.ContainsKey(key);
        }

        public T2 this[T1 key]
        {
            get
            {
                return forward[key];
            }
            set
            {
                forward[key] = value;
                reverse[value] = key;
            }
        }

        public T1 this[T2 key]
        {
            get
            {
                return reverse[key];
            }
            set
            {
                forward[value] = key;
                reverse[key] = value;
            }
        }
    }
}
