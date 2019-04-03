(import [pyhtml [*]])

(require [hy.contrib.walk [let]])
;; (require [build-templates.common [build-page]])

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

(defn form-required []
  (fieldset
    (legend "Search")

    (let [field-name "genome"]
      ((div :class_ "pure-control-group")
       ((label :for_ field-name)
        "Genome: ")

       ((select :id field-name :name field-name)
        ((option :value "dm6")
         "D. melanogaster"))))

    (let [field-name "sequences"]
      ((div :class_ "pure-control-group")
       ((label :for_ field-name)
        "Sequence(s): "
        ((span :class_ "pure-form-message")
         "One sequence per line"))

       (textarea
         :id field-name
         :name field-name
         :cols "30"
         :rows "2")))

    (let [field-name "submit"]
      ((div :class_ "pure-controls")
       ((button :id field-name :type "submit" :class_ "pure-button pure-button-primary")
        "Submit")))))

(defn build-query-page []
  (let [generated
        (build-page
          "Query"
          ((form :class_ "pure-form pure-form-aligned")
           (form-required)))]
    (+ "{% load static %}\n"
       (-> generated
           (.replace "fly.gif" "{% static \"search/fly.gif\" %}")
           (.replace "index.js" "{% static \"search/index.js\" %}")
           (.replace "pure-min.css" "{% static \"search/pure-min.css\" %}")
           (.replace "style.css" "{% static \"search/style.css\" %}")))))
