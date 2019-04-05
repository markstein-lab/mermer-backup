(import [pyhtml [*]])

(defmacro build-page [page-title &rest body-forms]
  "Emits the minimal code necessary for build a page.

`page-title' is used in the <title> tag."
  `(.render
     (html
       (head
         (meta :charset "utf-8")
         (title ~(.format "{} - Open Genomics : Genome Enhancer" page-title)))
       (body
         ((div :class_ "header")
          (img :src "fly.gif")
          (strong "merMER")
          (p "This program enables you to search entire genomes for clusters of DNA sequences"))
         ~@body-forms
         (script :src "index.js")
         (link :rel "stylesheet" :href "pure-min.css")
         (link :rel "stylesheet" :href "style.css")))))
