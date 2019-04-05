(import [os.path [dirname join]])

(import [search.query_page [build-query-page]])

(require [hy.contrib.walk [let]])

(setv *template-targets*
      [["templates/search/index.html" build-query-page]])

(defn expand-targets [targets]
  """Returns a list of pairs with the first element being the
absolute path for where the template should be written to, and
the second element being the function for generating the
template's contents.
"""
  (let [search-directory (dirname __file__)]
    (list (map (fn [target]
                 (let [path (first target)
                       build (second target)]
                   [(join search-directory path) build]))
               targets))))

(defn build-templates []
  """Construct Django template files and write them to the
relevant directory on disk.
"""
  (let [write-targets (expand-targets *template-targets*)]
    (for [target write-targets]
      (let [path (first target) build (second target)]
        (with [out (open path "w+")]
          (.write out (build)))))
    write-targets))
