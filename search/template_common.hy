(import [pyhtml [*]])

(require [hy.contrib.walk [let]])

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

(import [re [MULTILINE findall]])

(defn path-to-static [path-attribute]
  (let [url (.join "" (->> path-attribute
                           (drop-while (fn [c] (!= c "\"")))
                           (drop 1)
                           (drop-last 1)))]
    (+ (.join "" (take-while (fn [c] (!= c "\"")) path-attribute))
       "\"{% static \"search/"
       url
       "\" %}\"")))

(defn paths-to-static [html]
  (for [path (+ (findall r"src=\".*?\"" html MULTILINE)
                (findall r"href=\".*?\"" html MULTILINE))]
    (setv html (.replace html
                         path
                         (path-to-static path))))
  (+ "{% load static %}\n"
     html))
